# 04_coors_movement.R

RAW_FILE      <- "statcast_data/statcast_raw_2023_2026.rds"
ART_DIR       <- "app_artifacts"
OUT_DIR       <- "coors_out"
COORS         <- "COL"
PT            <- c("FF", "SI", "FC", "SL", "ST", "SV", "CU", "CH", "FS"
MIN_COORS     <- 10
MIN_ELSE      <- 50
MIN_GRADE     <- 100
GRADE_SEASONS <- 2025:2026
TARGET_SEASON <- 2026
TARGET_MIN    <- 400
RELABEL       <- TRUE
RETENTION     <- "by_type"
VELO_CREDIT   <- TRUE
N_BOOT        <- 200
SENS_CELLS    <- 120
SENS_K        <- 400
SEED          <- 42
THREADS       <- max(1L, parallel::detectCores() - 1L)
N_PROC        <- if (.Platform$OS.type == "windows") 1L else THREADS

pkgs <- c("data.table", "lightgbm", "ggplot2")
new  <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(new)) install.packages(new)
suppressPackageStartupMessages({ library(data.table); library(lightgbm); library(ggplot2) })
options(width = 220, scipen = 999, datatable.print.class = FALSE, datatable.showProgress = FALSE)
setDTthreads(THREADS); set.seed(SEED); t_start <- Sys.time()
dir.create(OUT_DIR, showWarnings = FALSE)
RPT  <- file.path(OUT_DIR, "coors_report.txt"); if (file.exists(RPT)) file.remove(RPT)
emit <- function(txt) { cat(txt, sep = "\n"); cat(txt, file = RPT, sep = "\n", append = TRUE) }
say  <- function(...) emit(paste0(...))
hdr  <- function(x) emit(c("", strrep("=", 110), x, strrep("=", 110)))
sub_ <- function(x) emit(c("", paste0("--- ", x, " ---")))
rnd  <- function(d, k = 3) d[, lapply(.SD, function(x) if (is.double(x)) round(x, k) else x)]
shw  <- function(x, k = 3) emit(capture.output(print(rnd(as.data.table(x), k), nrows = 500L, topn = 500L, row.names = FALSE)))
tick <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
wm   <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0; if (!any(ok)) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok]) }
wcor <- function(a, b, w) { ok <- is.finite(a) & is.finite(b) & is.finite(w) & w > 0
                            if (sum(ok) < 5) return(NA_real_); cov.wt(cbind(a[ok], b[ok]), wt = w[ok], cor = TRUE)$cor[1, 2] }
circ_mean <- function(a) { a <- a[is.finite(a)]; if (!length(a)) return(NA_real_)
                           (atan2(mean(sin(a * pi / 180)), mean(cos(a * pi / 180))) * 180 / pi) %% 360 }
adiff <- function(a, b) ((a - b + 180) %% 360) - 180
zs    <- function(x) { s <- sd(x, na.rm = TRUE); if (!is.finite(s) || s == 0) x * 0 else (x - mean(x, na.rm = TRUE)) / s }
slope0 <- function(yc, xe, w) { ok <- is.finite(yc) & is.finite(xe) & is.finite(w); sum(w[ok] * xe[ok] * yc[ok]) / sum(w[ok] * xe[ok]^2) }

# load
REQ <- c("game_type", "game_pk", "game_date", "game_year", "pitcher", "p_throws", "pitch_type", "release_speed", "release_spin_rate",
         "spin_axis", "release_extension", "release_pos_x", "release_pos_z", "pfx_x", "pfx_z", "plate_x", "plate_z",
         "vx0", "vy0", "vz0", "ax", "ay", "az", "home_team", "description")
OPT <- c("player_name", "arm_angle", "away_team", "inning_topbot", "bb_type", "launch_angle")
for (f in c(RAW_FILE, file.path(ART_DIR, c("grader.txt", "app_data.rds"))))
  if (!file.exists(f)) stop("Missing ", f, " (working directory ", getwd(), ")")
tick("loading ", RAW_FILE)
sc <- readRDS(RAW_FILE); setDT(sc)
miss <- setdiff(REQ, names(sc)); if (length(miss)) stop("Missing columns in ", RAW_FILE, ": ", paste(miss, collapse = ", "))
sc <- sc[game_type == "R", intersect(c(REQ, OPT), names(sc)), with = FALSE]; gc()
for (cc in setdiff(OPT, names(sc))) set(sc, j = cc, value = NA)
for (cc in c("release_speed", "release_spin_rate", "spin_axis", "release_extension", "release_pos_x", "release_pos_z", "pfx_x", "pfx_z",
             "plate_x", "plate_z", "vx0", "vy0", "vz0", "ax", "ay", "az", "arm_angle", "launch_angle"))
  if (!is.double(sc[[cc]])) set(sc, j = cc, value = suppressWarnings(as.numeric(sc[[cc]])))
sc[, `:=`(game_year = as.integer(game_year), game_date = as.IDate(game_date))]
sc[, pt := pitch_type][pitch_type %in% c("KC", "CS"), pt := "CU"]
n_raw <- nrow(sc)
sc <- sc[pt %in% PT & is.finite(release_speed) & is.finite(pfx_x) & is.finite(pfx_z) & is.finite(vx0) & is.finite(vy0) & is.finite(vz0) &
           is.finite(ax) & is.finite(ay) & is.finite(az) & ay > 0 & vy0 < 0 & is.finite(release_extension) & is.finite(release_pos_z) &
           is.finite(plate_x) & is.finite(plate_z)]
sc[, coors := home_team == COORS]
if (all(is.na(sc$player_name))) sc[, player_name := as.character(pitcher)]
sc[, rid := .I]

# per-pitch physics
GRAV <- 32.174; YP <- 17 / 12; K_SEA <- 5.383e-3
aero <- function(vx0, vy0, vz0, ax, ay, az, tp) {
  bx <- vx0 + ax * tp / 2; by <- vy0 + ay * tp / 2; bz <- vz0 + az * tp / 2; vb <- sqrt(bx^2 + by^2 + bz^2)
  aD <- -(ax * bx + ay * by + (az + GRAV) * bz) / vb
  list(vb = vb, ux = bx / vb, uy = by / vb, uz = bz / vb, aD = aD,
       tx = ax + aD * bx / vb, ty = ay + aD * by / vb, tz = az + GRAV + aD * bz / vb)
}
t_to <- function(vy, ay, dy) (-vy - sqrt(pmax(vy^2 - 2 * ay * dy, 0))) / ay      # time to travel dy (> 0 forward, < 0 back)
phys <- function(vx0, vy0, vz0, ax, ay, az, ext) {
  tp <- t_to(vy0, ay, 50 - YP); tr <- t_to(vy0, ay, 50 - (60.5 - ext))
  a  <- aero(vx0, vy0, vz0, ax, ay, az, tp); aT <- sqrt(a$tx^2 + a$ty^2 + a$tz^2)
  fx <- vx0 + ax * tp; fy <- vy0 + ay * tp; fz <- vz0 + az * tp
  list(flight = tp - tr, vplate = sqrt(fx^2 + fy^2 + fz^2) / 1.46667, kD = a$aD / a$vb^2, kT = aT / a$vb^2,
       vaa = -atan(fz / fy) * 180 / pi, cl = aT / (K_SEA * a$vb^2), vb = a$vb)
}
# The same pitch from the same release point to the same target, with lift scaled by bT and drag by bD
coors_sim <- function(vx0, vy0, vz0, ax, ay, az, ext, px, pz, bT, bD) {
  tp <- t_to(vy0, ay, 50 - YP); yr <- 60.5 - ext; tr <- t_to(vy0, ay, 50 - yr)
  x50 <- px - vx0 * tp - 0.5 * ax * tp^2; z50 <- pz - vz0 * tp - 0.5 * az * tp^2
  xr <- x50 + vx0 * tr + 0.5 * ax * tr^2; zr <- z50 + vz0 * tr + 0.5 * az * tr^2; vyr <- vy0 + ay * tr
  a <- aero(vx0, vy0, vz0, ax, ay, az, tp)
  ax2 <- bT * a$tx - bD * a$aD * a$ux; ay2 <- bT * a$ty - bD * a$aD * a$uy; az2 <- bT * a$tz - bD * a$aD * a$uz - GRAV
  t <- t_to(vyr, ay2, yr - YP)
  vx2 <- (px - xr - 0.5 * ax2 * t^2) / t; vz2 <- (pz - zr - 0.5 * az2 * t^2) / t
  fx <- vx2 + ax2 * t; fy <- vyr + ay2 * t; fz <- vz2 + az2 * t
  list(flight = t, vplate = sqrt(fx^2 + fy^2 + fz^2) / 1.46667, vaa = -atan(fz / fy) * 180 / pi)
}
tick("per-pitch physics on ", format(nrow(sc), big.mark = ","), " pitches")
ph <- sc[, phys(vx0, vy0, vz0, ax, ay, az, release_extension)]
sc[, c("flight", "vplate", "kD", "kT", "vaa", "cl", "vb") := ph]; rm(ph); gc()
sc[, sgn := fifelse(p_throws == "R", 1, -1)]
sc[, `:=`(hb_arm = -12 * pfx_x * sgn, ivb = 12 * pfx_z, axis_m = (180 + sgn * (spin_axis - 180)) %% 360)]   # arm side +, LHP mirrored
sc[, mov := sqrt(hb_arm^2 + ivb^2)]
# measured spin axis -> direction Magnus pushes the ball (from 12:00 toward arm side). Seam = the part across it.
sc[, th := (axis_m - 180) * pi / 180]
sc[, `:=`(m_mag = hb_arm * sin(th) + ivb * cos(th), m_seam = hb_arm * cos(th) - ivb * sin(th))]
# spin efficiency from the Magnus part only (Nathan's lift model, sea-level air), so seam movement doesn't inflate it
sc[, cl_m := cl * pmax(m_mag, 0) / pmax(mov, 0.1)]
sc[, eff_est := 78.92 * (0.4 * cl_m / (1 - 2.32 * cl_m)) * vb / release_spin_rate]
sc[!is.finite(eff_est) | cl_m >= 0.42 | eff_est > 1.3, eff_est := NA_real_]
sc[, c("cl", "cl_m", "vb") := NULL]
SWING <- c("swinging_strike", "swinging_strike_blocked", "foul", "foul_tip", "hit_into_play")
sc[, `:=`(swing = description %in% SWING, whiff = description %in% c("swinging_strike", "swinging_strike_blocked", "foul_tip"),
          take = description %in% c("ball", "blocked_ball", "called_strike"), cstr = description == "called_strike",
          bip = description == "hit_into_play")]
sc[, gb := bip & bb_type %in% "ground_ball"]
sc[bip == FALSE, launch_angle := NA_real_]

has_team <- !all(is.na(sc$inning_topbot)) && !all(is.na(sc$away_team))
sc[, p_team := if (has_team) fifelse(inning_topbot == "Top", home_team, away_team) else NA_character_]
PS <- sc[, .(pitches = .N, coors_share = mean(coors),
             team = if (all(is.na(p_team))) NA_character_ else names(which.max(table(p_team)))), by = .(pitcher, game_year)]
PS[, rockies := fifelse(is.na(team), coors_share >= 0.35, team == COORS)]
NAMES <- sc[, .(player_name = player_name[.N], p_throws = p_throws[.N]), by = pitcher]

# pitch classification
relab <- data.table(); relab_rate <- c(coors = 0, else_ = 0)
if (RELABEL) {
  okr <- sc[, is.finite(release_spin_rate) & is.finite(spin_axis)]
  cen <- sc[coors == FALSE & okr, .(cn = .N, cv = mean(release_speed), cs = mean(release_spin_rate), ca = circ_mean(axis_m)), by = .(pitcher, game_year, pt)][cn >= 20]
  sdv <- sc[coors == FALSE & okr, .(n = .N, v = sd(release_speed), s = sd(release_spin_rate), a = sd(adiff(axis_m, circ_mean(axis_m)))), by = .(pitcher, game_year, pt)][
    n >= 20, lapply(.SD, median, na.rm = TRUE), .SDcols = c("v", "s", "a")]
  best_d <- rep(Inf, nrow(sc)); best_pt <- sc$pt; own_d <- rep(Inf, nrow(sc))
  for (p in PT) {
    j <- cen[pt == p][sc, on = .(pitcher, game_year), .(cv = x.cv, cs = x.cs, ca = x.ca)]
    d2 <- ((sc$release_speed - j$cv) / sdv$v)^2 + ((sc$release_spin_rate - j$cs) / sdv$s)^2 + (adiff(sc$axis_m, j$ca) / sdv$a)^2
    d2[!is.finite(d2)] <- Inf
    own_d[sc$pt == p] <- d2[sc$pt == p]
    b <- d2 < best_d; best_d[b] <- d2[b]; best_pt[b] <- p
  }
  mv <- okr & best_pt != sc$pt & best_d < 9 & best_d < 0.5 * own_d
  relab <- data.table(coors = sc$coors[mv], from = sc$pt[mv], to = best_pt[mv])[, .(coors = sum(coors), elsewhere = sum(!coors)), by = .(from, to)]
  nfrom <- sc[, .(n_c = sum(coors), n_e = sum(!coors)), by = .(from = pt)]
  relab[nfrom, on = "from", `:=`(per1k_coors = 1000 * coors / i.n_c, per1k_else = 1000 * elsewhere / i.n_e)]
  setorder(relab, -coors)
  relab_rate <- c(coors = sum(mv & sc$coors) / sc[coors == TRUE, .N], else_ = sum(mv & !sc$coors) / sc[coors == FALSE, .N])
  sc[mv, pt := best_pt[mv]]
  rm(best_d, best_pt, own_d, d2, j, mv, okr, cen); gc()
}

# drag and lift relative to the same pitcher-pitch-season away from Coors
base <- sc[coors == FALSE, .(nE = .N, kD0 = median(kD, na.rm = TRUE), kT0 = median(kT, na.rm = TRUE)), by = .(pitcher, game_year, pt)][nE >= 20]
sc[base, on = .(pitcher, game_year, pt), `:=`(rD = kD / i.kD0, rT = kT / i.kT0)]
GM <- sc[is.finite(rD) & is.finite(rT) & rD > 0 & rT > 0,
         .(n = .N, D = median(rD), L = median(rT), park = home_team[1], game_year = game_year[1], month = month(game_date[1])), by = game_pk][n >= 50]

# same pitcher, pitch type, and season, at Coors and elsewhere
cell_stats <- function(d) d[, .(
  n = .N, velo = mean(release_speed), vplate = mean(vplate, na.rm = TRUE), flight = mean(flight, na.rm = TRUE),
  spin = mean(release_spin_rate, na.rm = TRUE), axis = circ_mean(axis_m), ivb = mean(ivb), hb = mean(hb_arm), mov = mean(mov),
  m_mag = mean(m_mag, na.rm = TRUE), m_seam = mean(m_seam, na.rm = TRUE),
  kT = mean(kT, na.rm = TRUE), kT_sd = sd(kT, na.rm = TRUE), kD = mean(kD, na.rm = TRUE), kD_sd = sd(kD, na.rm = TRUE),
  vaa = mean(vaa, na.rm = TRUE), relz = mean(release_pos_z), ext = mean(release_extension), arm = mean(arm_angle, na.rm = TRUE),
  eff = mean(eff_est, na.rm = TRUE), swings = sum(swing), whiffs = sum(whiff), takes = sum(take), cstr = sum(cstr),
  bips = sum(bip), gbs = sum(gb), la = mean(launch_angle, na.rm = TRUE)), by = .(pitcher, game_year, pt)]
tick("pairing Coors and non-Coors cells")
CE <- cell_stats(sc[coors == FALSE]); CC <- cell_stats(sc[coors == TRUE])
C  <- merge(CC, CE, by = c("pitcher", "game_year", "pt"), suffixes = c("_c", "_e"))[n_c >= MIN_COORS & n_e >= MIN_ELSE]
if (nrow(C) < 30) stop("Only ", nrow(C), " pitcher-season-pitch cells have ", MIN_COORS, "+ Coors pitches; lower MIN_COORS.")
C[PS, on = .(pitcher, game_year), `:=`(rockies = i.rockies, team = i.team)]
C[NAMES, on = "pitcher", `:=`(player_name = i.player_name, p_throws = i.p_throws)]
C[is.na(rockies), rockies := FALSE]
C[, w  := 1 / (1 / n_c + 1 / n_e)]                           
C[, `:=`(lT = log(kT_c / kT_e), lD = log(kD_c / kD_e))]
C[, `:=`(wT = 1 / (kT_sd_c^2 / (n_c * kT_c^2) + kT_sd_e^2 / (n_e * kT_e^2)),        # inverse variance of the log ratios
         wD = 1 / (kD_sd_c^2 / (n_c * kD_c^2) + kD_sd_e^2 / (n_e * kD_e^2)))]
C[, seam_share := abs(m_seam_e) / sqrt(m_mag_e^2 + m_seam_e^2)]

ret_fun <- function(D) {
  f <- function(d) list(
    cells = nrow(d), coors_pitches = sum(d$n_c),
    ivb = slope0(d$ivb_c, d$ivb_e, d$w), hb = slope0(d$hb_c, d$hb_e, d$w),
    magnus = slope0(d$m_mag_c, d$m_mag_e, d$w), seam = slope0(d$m_seam_c, d$m_seam_e, d$w),
    lift = exp(wm(d$lT, d$wT)), drag = exp(wm(d$lD, d$wD)),
    d_ivb = wm(d$ivb_c - d$ivb_e, d$w), d_hb = wm(d$hb_c - d$hb_e, d$w),
    d_velo = wm(d$velo_c - d$velo_e, d$w), d_vplate = wm(d$vplate_c - d$vplate_e, d$w),
    flight_pct = 100 * wm(d$flight_c / d$flight_e - 1, d$w), spin_pct = 100 * wm(d$spin_c / d$spin_e - 1, d$w),
    d_axis = wm(adiff(d$axis_c, d$axis_e), d$w), d_vaa = wm(d$vaa_c - d$vaa_e, d$w),
    d_relz = wm(d$relz_c - d$relz_e, d$w), d_ext = wm(d$ext_c - d$ext_e, d$w), d_arm = wm(d$arm_c - d$arm_e, d$w))
  rbind(D[, f(.SD), by = pt], D[, c(list(pt = "ALL"), f(.SD))])
}
boot_se <- function(D, FUN, keys, B = N_BOOT) {
  est <- FUN(D); ids <- unique(D$pitcher)
  reps <- rbindlist(lapply(seq_len(B), function(b) FUN(D[data.table(pitcher = sample(ids, length(ids), TRUE)), on = "pitcher", allow.cartesian = TRUE])), fill = TRUE)
  vars <- setdiff(names(est), keys)
  se <- reps[, lapply(.SD, function(x) sd(x[is.finite(x)])), by = keys, .SDcols = vars]
  list(est = est, se = se)
}
with_se <- function(b, vars, k = 2) {     # "est (se)" strings
  out <- copy(b$est)
  for (v in vars) { s <- b$se[out, on = names(out)[1], get(v)]; set(out, j = v, value = sprintf(paste0("%.", k, "f (%.", k, "f)"), out[[v]], s)) }
  out
}
tick("bootstrapping retention (", N_BOOT, " draws)")
R0 <- boot_se(C, ret_fun, "pt")

# retention used to build each pitch's Coors version
pool <- R0$est[pt == "ALL"]
shrink <- function(est, se, p) {        # empirical Bayes (DerSimonian-Laird) toward the precision-weighted mean
  ok <- is.finite(est) & is.finite(se) & se > 0
  if (sum(ok) < 3) return(fifelse(is.finite(est), est, p))
  w <- 1 / se[ok]^2; mu <- sum(w * est[ok]) / sum(w)
  tau2 <- max(0, (sum(w * (est[ok] - mu)^2) - (sum(ok) - 1)) / (sum(w) - sum(w^2) / sum(w)))
  out <- rep(mu, length(est)); out[ok] <- mu + tau2 / (tau2 + se[ok]^2) * (est[ok] - mu); out
}
RT <- data.table(pt = PT)[R0$est, on = "pt", `:=`(m = i.magnus, s = i.seam, l = i.lift, d = i.drag)]
RT[R0$se, on = "pt", `:=`(m_se = i.magnus, s_se = i.seam, l_se = i.lift, d_se = i.drag)]
if (RETENTION == "by_type") {
  RT[, `:=`(bM = shrink(m, m_se, pool$magnus), bS = shrink(s, s_se, pool$seam), bT = shrink(l, l_se, pool$lift), bD = shrink(d, d_se, pool$drag))]
} else RT[, `:=`(bM = pool$magnus, bS = pool$magnus, bT = pool$lift, bD = pool$drag)]
RT[!is.finite(bM), bM := pool$magnus][!is.finite(bS), bS := pool$seam][!is.finite(bT), bT := pool$lift][!is.finite(bD), bD := pool$drag]

# environment and movement
hdr("0. RUN INFO")
say("Run: ", format(Sys.time(), "%Y-%m-%d %H:%M"), " | ", R.version.string, " | data.table ", as.character(packageVersion("data.table")),
    " | lightgbm ", as.character(packageVersion("lightgbm")))
say("Regular-season rows: ", format(n_raw, big.mark = ","), " -> with a full 9-parameter fit and a modelled pitch type: ", format(nrow(sc), big.mark = ","),
    " | at Coors: ", format(sc[coors == TRUE, .N], big.mark = ","), " | seasons: ", paste(sort(unique(sc$game_year)), collapse = ", "))
say("Pitcher's team from inning_topbot/away_team: ", has_team, " | Rockies pitcher-seasons: ", PS[rockies == TRUE, .N])
say("Paired cells (", MIN_COORS, "+ Coors and ", MIN_ELSE, "+ elsewhere): ", nrow(C), " from ", uniqueN(C$pitcher), " pitchers | Rockies cells ",
    C[rockies == TRUE, .N], ", visitor cells ", C[rockies == FALSE, .N])
shw(C[, .(cells = .N, rockies = sum(rockies), coors_pitches = sum(n_c), median_coors_n = median(n_c)), by = pt][order(-cells)], 0)
if (RELABEL) {
  say("Re-labelled by velo, spin rate and spin axis: ", round(100 * relab_rate[["coors"]], 2), "% of Coors pitches, ",
      round(100 * relab_rate[["else_"]], 2), "% elsewhere. Higher at Coors, one way round = Savant labels drift there.")
  if (nrow(relab)) shw(head(relab, 12), 1)
}
emit("Reading guide: 'retention' = Coors value / value elsewhere for the same pitcher, pitch and season. With air density alone,")
emit("drag and lift both read the Coors density ratio (roughly 0.82), and movement in inches a little lower because the pitch arrives sooner.")

hdr("1. AIR AT EACH PARK, FROM DRAG (1.00 = the pitcher's own average away from Coors)")
pk <- GM[, .(games = .N, density_from_drag = wm(D, n), lift = wm(L, n)), by = park][order(density_from_drag)]
shw(pk, 3)
sub_("Coors by month and by season")
shw(GM[park == COORS, .(games = .N, density_from_drag = wm(D, n), lift = wm(L, n)), by = month][order(month)], 3)
shw(GM[park == COORS, .(games = .N, density_from_drag = wm(D, n), lift = wm(L, n)), by = game_year][order(game_year)], 3)
sub_("Lift vs drag: thin air alone moves them together, so lift/drag stays near 1")
pk[, lift_vs_drag := lift / density_from_drag]
oth <- pk[park != COORS]
say("Coors lift/drag ", round(pk[park == COORS, lift_vs_drag], 3), " | other parks: median ", round(median(oth$lift_vs_drag), 3),
    ", range ", round(min(oth$lift_vs_drag), 3), " to ", round(max(oth$lift_vs_drag), 3))
emit("Below that range = lift falls further than drag at Coors: the drag coefficient rising in thinner air (lower Reynolds")
emit("number), lift lost beyond density, or both. Section 3 splits lift into its Magnus and seam parts.")

hdr("2. WHAT CHANGES AT COORS: SAME PITCHER, SAME PITCH, SAME SEASON (estimate (bootstrap SE))")
emit("ivb/hb/magnus/seam: retention slopes, Coors on elsewhere. lift/drag: ratios of aerodynamic acceleration per v^2.")
emit("d_*: Coors minus elsewhere (velo/vplate mph; flight_pct and spin_pct in %; axis/arm degrees; relz/ext feet; vaa degrees, + = flatter).")
shw(with_se(R0, c("ivb", "hb", "magnus", "seam", "lift", "drag"), 3)[, .(pt, cells, coors_pitches, ivb, hb, magnus, seam, lift, drag)])
shw(with_se(R0, c("d_ivb", "d_hb", "d_velo", "d_vplate", "flight_pct", "spin_pct", "d_axis", "d_vaa"))[, .(pt, cells, d_ivb, d_hb, d_velo, d_vplate, flight_pct, spin_pct, d_axis, d_vaa)], 2)
shw(R0$est[, .(pt, d_relz, d_ext, d_arm)], 3)
sub_("Rockies vs visitors (all pitch types pooled): do pitchers who live there keep more?")
rv <- rbind(ret_fun(C[rockies == TRUE])[pt == "ALL"][, group := "Rockies"], ret_fun(C[rockies == FALSE])[pt == "ALL"][, group := "visitors"])
shw(rv[, .(group, cells, coors_pitches, ivb, hb, magnus, seam, lift, drag, d_velo, spin_pct, d_axis)], 3)

hdr("3. MAGNUS VS SEAM WAKE: WHICH PART OF THE MOVEMENT THIN AIR TAKES")
emit("magnus = movement along the direction the measured spin axis predicts; seam = movement across it (seam-shifted wake and other")
emit("non-Magnus forces). Equal retention = seam effects thin out exactly like Magnus. Seam inches are the elsewhere averages.")
s3 <- C[, .(cells = .N, magnus_in = wm(m_mag_e, w), seam_in = wm(m_seam_e, w), seam_share = wm(seam_share, w)), by = pt]
s3 <- merge(s3, with_se(R0, c("magnus", "seam"), 3)[, .(pt, magnus_ret = magnus, seam_ret = seam)], by = "pt")
shw(s3[order(-seam_share)], 3)
say("Retention used for each pitch's Coors version (", RETENTION, ", shrunk toward the pooled value):")
shw(RT[, .(pt, magnus = bM, seam = bS, lift = bT, drag = bD)], 3)

hdr("4. DOES ANY WAY OF MAKING MOVEMENT KEEP MORE OF IT? (cell-level lift retention)")
emit("Outcome: log(Coors lift / elsewhere lift) per cell. Traits are non-Coors values, standardised within pitch type, taken from")
emit("a different half of the pitches than the lift baseline. pct_per_sd = % change in the share of lift kept per 1 SD of the trait.")
emit("If air density were the whole story every row would be ~0: thin air takes the same share of lift however the spin is made.")
emit("Scale: 1% per SD moves retention from 0.80 to 0.81, about 0.15 in of four-seam ride; anything much under that is physics trivia.")
# split half: traits from one half of each cell's non-Coors pitches, the lift baseline from the other, so chance noise
# in the non-Coors sample can't sit on both sides of the regression
KEYS <- c("pitcher", "game_year", "pt")
sc[, half := rid %% 2L]
HA <- sc[coors == FALSE & half == 0L, .(nA = .N, kT_A = mean(kT, na.rm = TRUE), kT_sdA = sd(kT, na.rm = TRUE)), by = KEYS]
HB <- sc[coors == FALSE & half == 1L, .(spin_b = mean(release_spin_rate, na.rm = TRUE), eff_b = mean(eff_est, na.rm = TRUE),
                                        velo_b = mean(release_speed), arm_b = mean(arm_angle, na.rm = TRUE), relz_b = mean(release_pos_z),
                                        ext_b = mean(release_extension), mm = mean(m_mag, na.rm = TRUE), ss = mean(m_seam, na.rm = TRUE)), by = KEYS]
HB[, seam_b := abs(ss) / sqrt(mm^2 + ss^2)]
H <- merge(merge(C[, .(pitcher, game_year, pt, rockies, n_c, kT_c, kT_sd_c)], HA, by = KEYS), HB, by = KEYS)
H[, `:=`(lT = log(kT_c / kT_A), wT = 1 / (kT_sd_c^2 / (n_c * kT_c^2) + kT_sdA^2 / (nA * kT_A^2)))]
H <- H[is.finite(lT) & is.finite(wT) & wT > 0 & is.finite(eff_b) & is.finite(arm_b) & is.finite(seam_b)]
H[, `:=`(z_spin = zs(spin_b), z_eff = zs(eff_b), z_seam = zs(seam_b), z_velo = zs(velo_b), z_arm = zs(arm_b),
         z_relz = zs(relz_b), z_ext = zs(ext_b)), by = pt]
H[, wh := pmin(wT, quantile(wT, .99))]
TERMS <- c("z_spin", "z_eff", "z_seam", "z_velo", "z_arm", "z_relz", "z_ext", "rockiesTRUE")
het_form <- as.formula(paste("lT ~", paste(c(if (uniqueN(H$pt) > 1) "pt", if (uniqueN(H$game_year) > 1) "factor(game_year)",
                                              "z_spin + z_eff + z_seam + z_velo + z_arm + z_relz + z_ext", if (uniqueN(H$rockies) > 1) "rockies"), collapse = " + ")))
het_fun <- function(D) { cf <- coef(lm(het_form, data = D, weights = wh)); data.table(term = TERMS, b = unname(cf[TERMS])) }
if (nrow(H) >= 60) {
  hb <- boot_se(H, het_fun, "term")
  ht <- merge(hb$est, hb$se[, .(term, se = b)], by = "term")[match(TERMS, term)]
  ht[, `:=`(pct_per_sd = 100 * (exp(b) - 1), ci_lo = 100 * (exp(b - 1.96 * se) - 1), ci_hi = 100 * (exp(b + 1.96 * se) - 1))]
  ht[, term := c("spin rate", "spin efficiency (est.)", "seam share", "velocity", "arm angle", "release height", "extension", "Rockies pitcher")]
  say("Cells: ", nrow(H), " | pooled lift retention ", round(exp(wm(H$lT, H$wh)), 3))
  shw(ht[, .(trait = term, pct_per_sd, ci_lo, ci_hi)], 2)
} else say("Too few cells (", nrow(H), ") for the trait model.")

# Pitch Shaper grader on sea-level and Coors versions of every pitch
tick("building Coors versions of ", GRADE_SEASONS[1], "-", max(GRADE_SEASONS), " pitches")
E <- sc[coors == FALSE & game_year %in% GRADE_SEASONS]
E[RT, on = "pt", `:=`(bM = i.bM, bS = i.bS, bT = i.bT, bD = i.bD)]
E[, `:=`(hb_c  = fifelse(is.finite(m_mag), bM * m_mag * sin(th) + bS * m_seam * cos(th), bM * hb_arm),
         ivb_c = fifelse(is.finite(m_mag), bM * m_mag * cos(th) - bS * m_seam * sin(th), bM * ivb))]
s1 <- E[, coors_sim(vx0, vy0, vz0, ax, ay, az, release_extension, plate_x, plate_z, bT, bD)]
s0 <- E[, coors_sim(vx0, vy0, vz0, ax, ay, az, release_extension, plate_x, plate_z, 1, 1)]
E[, `:=`(flight_cs = s1$flight, vplate_cs = s1$vplate, vaa_cs = s1$vaa, vaa_bs = s0$vaa, velo_eq = release_speed * s0$flight / s1$flight)]
rm(s0, s1)
vfit <- E[is.finite(vaa_bs), as.list(coef(lm(vaa_bs ~ plate_z))), by = pt]; setnames(vfit, c("pt", "v0", "v1"))
E[vfit, on = "pt", `:=`(vaa_x = vaa_bs - (i.v0 + i.v1 * plate_z), vaa_x_c = vaa_cs - (i.v0 + i.v1 * plate_z))]   # flatness vs league at that height

GC <- E[, .(n = .N, release_speed = mean(release_speed), hb_arm = mean(hb_arm), ivb = mean(ivb), release_extension = mean(release_extension),
            arm_angle = mean(arm_angle, na.rm = TRUE), release_pos_z = mean(release_pos_z), hb_c = mean(hb_c), ivb_c = mean(ivb_c),
            velo_eq = mean(velo_eq, na.rm = TRUE), vplate = mean(vplate, na.rm = TRUE), vplate_c = mean(vplate_cs, na.rm = TRUE),
            vaa_x = mean(vaa_x, na.rm = TRUE), vaa_x_c = mean(vaa_x_c, na.rm = TRUE), spin = mean(release_spin_rate, na.rm = TRUE),
            eff = mean(eff_est, na.rm = TRUE), m_mag = mean(m_mag, na.rm = TRUE), m_seam = mean(m_seam, na.rm = TRUE), mov = mean(mov)),
        by = .(pitcher, game_year, pt)]
rm(E); gc()
GC[NAMES, on = "pitcher", `:=`(player_name = i.player_name, p_throws = i.p_throws)]
GC <- GC[n >= MIN_GRADE & is.finite(arm_angle) & is.finite(release_extension) & is.finite(release_pos_z)]
GC[, `:=`(mov_c = sqrt(hb_c^2 + ivb_c^2), seam_share = abs(m_seam) / sqrt(m_mag^2 + m_seam^2))]
fb <- GC[pt %in% c("FF", "SI", "FC")][order(-n), .(fb_type = pt[1], fb_velo = release_speed[1], fb_hb = hb_arm[1], fb_ivb = ivb[1],
                                                     fb_velo_eq = velo_eq[1], fb_hb_c = hb_c[1], fb_ivb_c = ivb_c[1]), by = .(pitcher, game_year)]
GC[fb, on = .(pitcher, game_year), `:=`(fb_type = i.fb_type, fb_velo = i.fb_velo, fb_hb = i.fb_hb, fb_ivb = i.fb_ivb,
                                         fb_velo_eq = i.fb_velo_eq, fb_hb_c = i.fb_hb_c, fb_ivb_c = i.fb_ivb_c)]
GC[, is_fb := as.numeric(!is.na(fb_type) & pt == fb_type)]
GC[, is_lhp := as.numeric(p_throws == "L")]

A  <- readRDS(file.path(ART_DIR, "app_data.rds")); GB <- lgb.load(file.path(ART_DIR, "grader.txt"))
GF <- A$grader_feats; STUFF <- A$stuff; BM <- A$bat_means
grade_block <- function(shapes, pool, same_hand) {
  K <- nrow(pool$ctx); chunk <- max(1L, min(floor(150000 / K), ceiling(nrow(shapes) / N_PROC)))
  chunks <- split(seq_len(nrow(shapes)), ceiling(seq_len(nrow(shapes)) / chunk))
  one <- function(ix) {
    s <- rep(ix, each = K); r <- rep(seq_len(K), length(ix))
    X <- matrix(NA_real_, length(s), length(GF), dimnames = list(NULL, GF))
    for (v in STUFF) if (v %in% GF) X[, v] <- shapes[[v]][s] + pool$dev[r, v]
    fbm <- shapes$is_fb[s] == 1
    dif <- function(v, ref) ifelse(fbm, 0, shapes[[v]][s] + pool$dev[r, v] - shapes[[ref]][s])
    if ("velo_diff" %in% GF) X[, "velo_diff"] <- dif("release_speed", "fb_velo")
    if ("hb_diff" %in% GF)   X[, "hb_diff"]   <- dif("hb_arm", "fb_hb")
    if ("ivb_diff" %in% GF)  X[, "ivb_diff"]  <- dif("ivb", "fb_ivb")
    if ("is_fb" %in% GF)     X[, "is_fb"]     <- shapes$is_fb[s]
    if ("is_lhp" %in% GF)    X[, "is_lhp"]    <- shapes$is_lhp[s]
    if ("same_hand" %in% GF) X[, "same_hand"] <- same_hand
    for (v in intersect(c("px_bat", "pz_norm", "balls", "strikes"), GF)) X[, v] <- pool$ctx[r, v]
    for (v in intersect(names(BM), GF)) X[, v] <- BM[[v]]
    100 * colMeans(matrix(if (N_PROC > 1L) predict(GB, X, params = list(num_threads = 1L)) else predict(GB, X), nrow = K))
  }
  res <- if (N_PROC > 1L && length(chunks) > 1L) parallel::mclapply(chunks, one, mc.cores = N_PROC) else lapply(chunks, one)
  bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
  if (any(bad)) res[bad] <- lapply(chunks[bad], one)
  unlist(res, use.names = FALSE)
}
grade_xrv <- function(S, pools = A$pool) {
  out <- matrix(NA_real_, nrow(S), 2)
  for (h in 0:1) for (p in unique(S$pt)) {
    i <- which(S$pt == p); kk <- paste(p, h); if (is.null(pools[[kk]])) next
    out[i, h + 1L] <- grade_block(S[i], pools[[kk]], h) + unname(A$offset[paste(p, h, S$is_lhp[i])])
  }
  rowMeans(out)
}
to_plus <- function(x) 100 - 10 * (x - A$scale[["mu"]]) / A$scale[["sd"]]
shape <- function(version, d = GC) {
  s <- d[, .(pt, is_lhp, release_speed, hb_arm, ivb, release_extension, arm_angle, release_pos_z, fb_velo, fb_hb, fb_ivb, is_fb)]
  if (version != "sea") s[, `:=`(hb_arm = d$hb_c, ivb = d$ivb_c, fb_hb = d$fb_hb_c, fb_ivb = d$fb_ivb_c)]
  if (version == "full") s[, `:=`(release_speed = d$velo_eq, fb_velo = d$fb_velo_eq)]
  s
}
CKPT <- file.path(OUT_DIR, "graded_cells.rds")
ck <- if (file.exists(CKPT)) tryCatch(readRDS(CKPT), error = function(e) NULL) else NULL
same <- !is.null(ck) && nrow(ck) == nrow(GC) && isTRUE(VELO_CREDIT == attr(ck, "velo_credit")) &&
  all(ck$pitcher == GC$pitcher & ck$game_year == GC$game_year & ck$pt == GC$pt) &&
  isTRUE(all.equal(c(ck$hb_c, ck$ivb_c, ck$velo_eq, ck$fb_hb_c), c(GC$hb_c, GC$ivb_c, GC$velo_eq, GC$fb_hb_c)))
if (same) {
  tick("reusing grades from ", CKPT); GC[, c("g_sea", "g_move", "g_full") := ck[, .(g_sea, g_move, g_full)]]
} else {
  tick("grading ", nrow(GC), " pitcher-season-pitch shapes (sea level and Coors) on ", N_PROC, " processes")
  GC[, g_sea  := to_plus(grade_xrv(shape("sea")))]
  GC[, g_move := to_plus(grade_xrv(shape("move")))]
  GC[, g_full := if (VELO_CREDIT) to_plus(grade_xrv(shape("full"))) else g_move]
  ck <- GC[, .(pitcher, game_year, pt, hb_c, ivb_c, velo_eq, fb_hb_c, g_sea, g_move, g_full)]; setattr(ck, "velo_credit", VELO_CREDIT)
  saveRDS(ck, CKPT)
}
rm(ck)
GC[, `:=`(loss = g_sea - g_full, loss_move = g_sea - g_move)]

hdr(paste0("5. VALUE AT COORS: PITCH SHAPER GRADES OF EACH PITCH'S COORS VERSION (", paste(range(GRADE_SEASONS), collapse = "-"), ", ", MIN_GRADE, "+ pitches elsewhere)"))
emit("Same release, velocity, location mix and batters; movement replaced by its Coors version. sea = grade elsewhere,")
emit(paste0("coors_move = movement change only, coors = movement plus the velocity credit for a shorter flight (", if (VELO_CREDIT) "on" else "off", ")."))
emit("loss = sea - coors in grade points (10 = 1 SD of all pitches). Pitch-weighted.")
t5 <- GC[, .(cells = .N, sea = wm(g_sea, n), coors_move = wm(g_move, n), coors = wm(g_full, n), loss = wm(loss, n),
             loss_sd = sd(loss), pct_losing_5plus = 100 * wm(as.numeric(loss >= 5), n),
             d_ivb = wm(ivb_c - ivb, n), d_hb = wm(hb_c - hb_arm, n), velo_credit = wm(velo_eq - release_speed, n),
             d_flatness = wm(vaa_x_c - vaa_x, n)), by = pt][order(loss)]
shw(t5, 2)
say("All pitches: sea ", round(GC[, wm(g_sea, n)], 1), " | Coors ", round(GC[, wm(g_full, n)], 1), " | loss ", round(GC[, wm(loss, n)], 2),
    " grade points | movement alone ", round(GC[, wm(loss_move, n)], 2), ", velocity credit gives back ", round(GC[, wm(g_full - g_move, n)], 2))

hdr("6. SAME SHAPE, DIFFERENT RECIPE: WHAT EACH TRAIT IS WORTH AT SEA LEVEL VS AT COORS")
emit("Grade points from +1 SD (within pitch type) of one trait with everything else held, graded elsewhere and at Coors.")
emit("movement = +1 SD of total movement in the pitch's own direction, thinned at Coors like the rest of its movement.")
emit("Spin rate and spin efficiency don't enter the grade directly; they act only through movement, and section 4 shows")
emit(paste0("whether they change how much of it survives. Up to ", SENS_CELLS, " cells per pitch type; kept = Coors worth / sea-level worth."))
SV_T <- c(velocity = "release_speed", extension = "release_extension", `release height` = "release_pos_z", `arm angle` = "arm_angle", movement = "mov")
sdp  <- GC[, lapply(.SD, sd, na.rm = TRUE), by = pt, .SDcols = unname(SV_T)]
smp  <- merge(GC[is.finite(velo_eq), .SD[sample(.N, min(.N, SENS_CELLS))], by = pt],
              setnames(copy(sdp), c("pt", paste0("sd_", unname(SV_T)))), by = "pt")
POOL_S <- lapply(A$pool, function(pl) { k <- sort(sample(nrow(pl$ctx), min(SENS_K, nrow(pl$ctx)))); list(ctx = pl$ctx[k, , drop = FALSE], dev = pl$dev[k, , drop = FALSE]) })
tick("trait worth at sea level vs Coors on ", nrow(smp), " cells")
sens <- tryCatch({
  b_sea <- shape("sea", smp); b_co <- shape("full", smp)
  g0s <- to_plus(grade_xrv(b_sea, POOL_S)); g0c <- to_plus(grade_xrv(b_co, POOL_S))
  rbindlist(lapply(names(SV_T), function(k) {
    v <- SV_T[[k]]; s1 <- copy(b_sea); c1 <- copy(b_co); d <- smp[[paste0("sd_", v)]]
    if (v == "mov") { f <- 1 + d / pmax(smp$mov, 1); s1[, `:=`(hb_arm = hb_arm * f, ivb = ivb * f)]; c1[, `:=`(hb_arm = hb_arm * f, ivb = ivb * f)] }
    else { set(s1, j = v, value = s1[[v]] + d); set(c1, j = v, value = c1[[v]] + d) }
    data.table(pt = smp$pt, trait = k, sea = to_plus(grade_xrv(s1, POOL_S)) - g0s, coors = to_plus(grade_xrv(c1, POOL_S)) - g0c)
  }))
}, error = function(e) { say("Section 6 failed: ", conditionMessage(e)); NULL })
if (!is.null(sens)) {
  st <- sens[, .(sea = mean(sea, na.rm = TRUE), coors = mean(coors, na.rm = TRUE)), by = .(pt, trait)]
  st[, cell := sprintf("%+.2f / %+.2f", sea, coors)]
  shw(dcast(st, pt ~ factor(trait, levels = names(SV_T)), value.var = "cell"))
  emit("cells: worth at sea level / worth at Coors, grade points per +1 SD. The SDs:")
  shw(setnames(copy(sdp), c("pt", names(SV_T))), 2)
  sub_("All pitch types: worth per +1 SD")
  shw(sens[, .(sea = mean(sea, na.rm = TRUE), coors = mean(coors, na.rm = TRUE)), by = trait][, kept := coors / sea], 2)
}
sub_("Four-seam flatness: where it comes from vs what Coors leaves (flatness = VAA above league at that plate height, degrees)")
ff <- GC[pt == "FF"]
if (nrow(ff) >= 60) {
  ff[, `:=`(height = cut(release_pos_z, quantile(release_pos_z, 0:3 / 3), labels = c("low release", "mid release", "high release"), include.lowest = TRUE),
            ride = cut(ivb, quantile(ivb, 0:3 / 3), labels = c("low IVB", "mid IVB", "high IVB"), include.lowest = TRUE))]
  shw(dcast(ff[, .(v = sprintf("%+.2f -> %+.2f (%.1f)", wm(vaa_x, n), wm(vaa_x_c, n), wm(loss, n))), by = .(height, ride)], height ~ ride, value.var = "v"))
  emit("cells: flatness elsewhere -> at Coors (grade loss)")
  fm <- lm(I(vaa_x_c - vaa_x) ~ ivb + release_pos_z + release_extension + release_speed, data = ff, weights = n)
  say("Change in flatness at Coors per unit: IVB ", round(coef(fm)[["ivb"]], 3), " deg/in | release height ", round(coef(fm)[["release_pos_z"]], 3),
      " deg/ft | extension ", round(coef(fm)[["release_extension"]], 3), " deg/ft | velo ", round(coef(fm)[["release_speed"]], 3), " deg/mph")
}

hdr("7. CHECK AGAINST WHAT ACTUALLY HAPPENS AT COORS (paired cells, all seasons)")
emit("Coors minus elsewhere, same pitcher-pitch-season. whiff = whiffs/swings, cstr = called strikes/takes, gb = ground balls/balls in play,")
emit("la = launch angle on balls in play. Hitters differ too (Rockies hitters at home, visitors used to more movement), so read as a check, not a test.")
out_fun <- function(D) {
  f <- function(d) {
    ws <- 1 / (1 / pmax(d$swings_c, 1) + 1 / pmax(d$swings_e, 1)); wt <- 1 / (1 / pmax(d$takes_c, 1) + 1 / pmax(d$takes_e, 1))
    wb <- 1 / (1 / pmax(d$bips_c, 1) + 1 / pmax(d$bips_e, 1))
    list(cells = nrow(d), d_whiff = 100 * wm(d$whiffs_c / d$swings_c - d$whiffs_e / d$swings_e, ws),
         d_cstr = 100 * wm(d$cstr_c / d$takes_c - d$cstr_e / d$takes_e, wt), d_gb = 100 * wm(d$gbs_c / d$bips_c - d$gbs_e / d$bips_e, wb),
         d_la = wm(d$la_c - d$la_e, wb))
  }
  rbind(D[, f(.SD), by = pt], D[, c(list(pt = "ALL"), f(.SD))])
}
OB <- boot_se(C, out_fun, "pt")
o7 <- merge(with_se(OB, c("d_whiff", "d_cstr", "d_gb", "d_la"), 1), t5[, .(pt, predicted_loss = loss, d_ivb_model = d_ivb)], by = "pt", all.x = TRUE)
shw(o7[order(predicted_loss)], 2)
pt_rows <- merge(OB$est[pt != "ALL"], t5[, .(pt, predicted_loss = loss)], by = "pt")
say("Across pitch types (weighted by cells): r(predicted grade loss, whiff change) = ", round(pt_rows[, wcor(predicted_loss, d_whiff, cells)], 2),
    "  (negative = types predicted to lose more do lose more whiffs)")
cl <- merge(C, GC[, .(pitcher, game_year, pt, loss)], by = c("pitcher", "game_year", "pt"))
if (nrow(cl) >= 30)
  say("Cell level (", nrow(cl), " graded cells): r(predicted loss, whiff change) = ",
      round(cl[, wcor(loss, whiffs_c / pmax(swings_c, 1) - whiffs_e / pmax(swings_e, 1), 1 / (1 / pmax(swings_c, 1) + 1 / pmax(swings_e, 1)))], 3))

# pitcher archetypes and targets
GC[, use := n / sum(n), by = .(pitcher, game_year)]
GRP <- c(FF = "second fastball", SI = "second fastball", FC = "cutter/gyro slider", SL = "cutter/gyro slider", ST = "sweeper/curve", SV = "sweeper/curve",
         CU = "sweeper/curve", CH = "changeup/splitter", FS = "changeup/splitter")
PP <- GC[, .(pitches = sum(n), types = .N, sea = sum(use * g_sea), coors = sum(use * g_full), coors_move = sum(use * g_move),
             fb_type = fb_type[1], fb_velo = fb_velo[1], fb_ivb = fb_ivb[1], fb_hb = fb_hb[1],
             fb_flat = if (any(is_fb == 1)) vaa_x[is_fb == 1][1] else NA_real_,
             arm_angle = sum(use * arm_angle), release_pos_z = sum(use * release_pos_z), release_extension = sum(use * release_extension),
             use_sweep_curve = sum(use[pt %in% c("ST", "SV", "CU")]), use_gyro = sum(use[pt == "SL" | (pt == "FC" & is_fb == 0)]),
             use_off = sum(use[pt %in% c("CH", "FS")]), mov = sum(use * mov), seam = sum(use * seam_share, na.rm = TRUE),
             main_sec = { o <- which(is_fb == 0); if (length(o)) GRP[[pt[o][which.max(n[o])]]] else "none" }),
         by = .(pitcher, player_name, game_year, p_throws)]
PP[, loss := sea - coors]
PP[PS, on = .(pitcher, game_year), `:=`(team = i.team, rockies = i.rockies)]
PP[, fb_family := fcase(is.na(fb_type), "no fastball", fb_type == "FF" & fb_ivb >= 17, "four-seam, 17+ in IVB", fb_type == "FF", "four-seam, under 17 in",
                        fb_type == "SI", "sinker", default = "cutter")]
PP[, slot := cut(arm_angle, c(-Inf, 25, 40, 52, Inf), labels = c("low/side (<25)", "3/4 low (25-40)", "3/4 high (40-52)", "over the top (52+)"))]
PQ <- PP[pitches >= TARGET_MIN]

hdr(paste0("8. PITCHER ARCHETYPES (", paste(range(GRADE_SEASONS), collapse = "-"), ", ", TARGET_MIN, "+ pitches elsewhere; arsenal grade = usage-weighted)"))
arch <- function(by) PQ[, .(pitchers = .N, sea = wm(sea, pitches), coors = wm(coors, pitches), loss = wm(loss, pitches),
                            loss_sd = sd(loss)), by = by][order(loss)]
for (b in c("fb_family", "slot", "main_sec")) { sub_(paste("by", b)); shw(arch(b), 2) }
sub_("primary fastball x main secondary (8+ pitcher-seasons)")
shw(arch(c("fb_family", "main_sec"))[pitchers >= 8], 2)
sub_("Which arsenal traits predict the loss? Grade points per 1 SD (bootstrap SE, resampling pitchers)")
PT_V <- c("fb_velo", "fb_ivb", "fb_hb", "arm_angle", "release_pos_z", "release_extension", "use_sweep_curve", "use_gyro", "use_off", "mov", "seam")
pm_fun <- function(D) { D <- copy(D); for (v in PT_V) set(D, j = v, value = zs(D[[v]]))
  cf <- coef(lm(as.formula(paste("loss ~", paste(PT_V, collapse = " + "))), data = D, weights = pmin(pitches, 3000)))
  data.table(term = PT_V, b = unname(cf[PT_V])) }
PQ2 <- PQ[complete.cases(PQ[, PT_V, with = FALSE])]
if (nrow(PQ2) >= 50) {
  pb <- boot_se(PQ2, pm_fun, "term")
  shw(merge(pb$est, pb$se[, .(term, se = b)], by = "term")[match(PT_V, term)][, .(trait = term, grade_pts_per_sd = b, se)], 2)
  emit("mov = usage-weighted total movement (in), seam = usage-weighted seam share, use_* = usage shares; fb_hb is arm-side run.")
}

hdr(paste0("9. TARGETS (", TARGET_SEASON, ", ", TARGET_MIN, "+ pitches elsewhere)"))
TT <- PQ[game_year == TARGET_SEASON]
if (nrow(TT) >= 20) {
  ef <- lm(coors ~ sea, data = TT, weights = pitches)
  TT[, `:=`(edge = coors - predict(ef, TT), rockies_value = (sea + coors) / 2)]
  emit("edge = Coors grade minus what pitchers of the same sea-level grade keep at Coors: how much better than his market")
  emit("(sea-level) value a pitcher plays at Coors. rockies_value = half home, half road.")
  cols <- c("player_name", "team", "p_throws", "pitches", "fb_type", "fb_velo", "fb_ivb", "arm_angle", "main_sec", "sea", "coors", "loss", "edge", "rockies_value")
  sub_("Non-Rockies, league average or better at sea level, ranked by edge")
  shw(head(TT[rockies %in% FALSE & sea >= 100][order(-edge), cols, with = FALSE], 25), 1)
  sub_("Non-Rockies, ranked by Rockies value")
  shw(head(TT[rockies %in% FALSE][order(-rockies_value), cols, with = FALSE], 25), 1)
  sub_("Current Rockies")
  shw(TT[rockies %in% TRUE][order(-rockies_value), cols, with = FALSE], 1)
} else say("Fewer than 20 qualifying pitchers in ", TARGET_SEASON, ".")

hdr("10. HOW PITCHERS ADAPT: PITCH MIX AT COORS (pitcher-seasons with 60+ Coors and 200+ other pitches)")
U <- sc[, .N, by = .(pitcher, game_year, pt, coors)]
U <- dcast(U, pitcher + game_year + pt ~ coors, value.var = "N", fill = 0); setnames(U, c("FALSE", "TRUE"), c("n_e", "n_c"))
U[, `:=`(tot_e = sum(n_e), tot_c = sum(n_c)), by = .(pitcher, game_year)]
U <- U[tot_c >= 60 & tot_e >= 200]
U[PS, on = .(pitcher, game_year), rockies := i.rockies]
shw(dcast(U[, .(d_share = 100 * wm(n_c / tot_c - n_e / tot_e, tot_c)), by = .(pt, group = fifelse(rockies %in% TRUE, "Rockies", "visitors"))],
          pt ~ group, value.var = "d_share"), 2)
emit("Percentage-point change in usage at Coors (pitch-type labels after re-labelling).")

# visualization
fwrite(pk, file.path(OUT_DIR, "parks.csv"))
fwrite(C, file.path(OUT_DIR, "cells_paired.csv"))
fwrite(GC[, .(pitcher, player_name, game_year, pt, p_throws, n, release_speed, hb_arm, ivb, hb_c, ivb_c, velo_eq, spin, eff, seam_share,
              arm_angle, release_pos_z, release_extension, vaa_x, vaa_x_c, g_sea, g_move, g_full, loss)], file.path(OUT_DIR, "pitch_grades.csv"))
fwrite(PP, file.path(OUT_DIR, "pitchers.csv"))
if (exists("TT") && "edge" %in% names(TT)) fwrite(TT[order(-edge)], file.path(OUT_DIR, paste0("targets_", TARGET_SEASON, ".csv")))

tick("plots")
rl <- melt(R0$est[pt != "ALL", .(pt, ivb, hb, magnus, seam, lift, drag)], id.vars = "pt", value.name = "est")
rl <- merge(rl, melt(R0$se[pt != "ALL", .(pt, ivb, hb, magnus, seam, lift, drag)], id.vars = "pt", value.name = "se"), by = c("pt", "variable"))
p1 <- ggplot(rl[is.finite(est)], aes(x = pt, y = est, ymin = est - 1.96 * se, ymax = est + 1.96 * se)) +
  geom_hline(yintercept = pool$drag, linetype = 2, colour = "grey50") + geom_pointrange(size = 0.3) + facet_wrap(~ variable, nrow = 2) +
  coord_cartesian(ylim = c(0.4, 1.2)) + labs(x = NULL, y = "Coors / elsewhere", title = "Retention at Coors by pitch type (dashed: drag ratio)") + theme_minimal(11)
ggsave(file.path(OUT_DIR, "01_retention_by_type.png"), p1, width = 11, height = 6, dpi = 110)
p2 <- ggplot(C, aes(x = ivb_e, y = ivb_c, colour = pt, size = w)) + geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_abline(slope = pool$ivb, intercept = 0, colour = "grey30") + geom_point(alpha = 0.5) + scale_size_area(max_size = 3, guide = "none") +
  labs(x = "IVB elsewhere (in)", y = "IVB at Coors (in)", colour = NULL, title = "Same pitcher, pitch and season: IVB at Coors vs elsewhere") + theme_minimal(11)
ggsave(file.path(OUT_DIR, "02_ivb_coors_vs_elsewhere.png"), p2, width = 8, height = 7, dpi = 110)
p3 <- ggplot(GC, aes(x = reorder(pt, loss), y = loss)) + geom_hline(yintercept = 0, colour = "grey60") + geom_boxplot(outlier.size = 0.4) +
  labs(x = NULL, y = "Grade points lost at Coors", title = "Pitch Shaper grade lost at Coors, by pitch type") + theme_minimal(11)
ggsave(file.path(OUT_DIR, "03_grade_loss_by_type.png"), p3, width = 8, height = 5, dpi = 110)
if (exists("TT") && "edge" %in% names(TT)) {
  lab <- rbind(head(TT[rockies %in% FALSE & sea >= 100][order(-edge)], 12), TT[rockies %in% TRUE])
  p4 <- ggplot(TT, aes(x = sea, y = coors)) + geom_abline(slope = coef(ef)[2], intercept = coef(ef)[1], colour = "grey50") +
    geom_point(aes(colour = rockies %in% TRUE), alpha = 0.6) + scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "#33006F"), labels = c("other", "Rockies"), name = NULL) +
    geom_text(data = lab, aes(label = player_name), size = 2.6, vjust = -0.7, check_overlap = TRUE) +
    labs(x = "Arsenal grade elsewhere", y = "Arsenal grade at Coors", title = paste0(TARGET_SEASON, " arsenals: above the line keeps more at Coors")) + theme_minimal(11)
  ggsave(file.path(OUT_DIR, "04_pitchers_sea_vs_coors.png"), p4, width = 9, height = 7, dpi = 110)
}

hdr("11. DONE")
say("Files in ", OUT_DIR, ": ", paste(list.files(OUT_DIR), collapse = ", "))
say("Elapsed (min): ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1))
