#!/usr/bin/env Rscript
# CEM (fia-cem-maine engine) per plot archive -> gate_components.csv and trajectory.csv. cem_v3 converter (track E F3).
# Tool side quantities and the expansion are a port of the registered scorer cem_gate_score_s41d8.R
# (version A point estimate, mean over sims; the bootstrap interval is left to gate_score). The only
# change: the plot frame of the EXPVOL evaluation is restricted to the AOI's (STATECD, COUNTYCD) when
# --counties is given (attribute selection; no coordinates are read).
# Usage: Rscript convert_cem.R --archive <per_plot_projections.rds> --state ME --fia_dir <dir> --work <dir>
#        --donors ME,NH,VT,NY,MA,CT,RI --eval_year 2024 [--counties 23019,23021] --out <dir>
#        [--aoi_id x --arm x --run_id x --engine_version x --config_md5 x]
#        [--scen_archives "bau:rcp45=<rds>;no_harvest:rcp45=<rds>"]  (trajectory rows only)
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); o <- list(counties = "", eval_year = "2024", aoi_id = "", arm = "", run_id = "",
                                  engine_version = "", config_md5 = "")
for (i in seq(1, length(a), 2)) o[[sub("^--", "", a[i])]] <- a[i + 1]
ST <- o$state; LB_TG <- 4.53592e-10
SC <- c(ME = 23, NH = 33, VT = 50, MI = 26, WI = 55, MN = 27, NY = 36, MA = 25, CT = 9, RI = 44, IA = 19,
        IL = 17, OH = 39, IN = 18, GA = 13)
sc <- SC[[ST]]; f <- function(st, tb) file.path(o$fia_dir, sprintf("%s_%s.csv", st, tb))
cty <- if (nzchar(o$counties)) as.integer(substr(strsplit(o$counties, ",")[[1]], 3, 5)) else integer()

tool <- function(path) {
  d <- as.data.table(readRDS(path))[STATECD == sc & scenario == "BAU"]
  if (length(cty)) d <- d[COUNTYCD %in% cty]
  d[, T2_PLT_CN := as.character(T2_PLT_CN)]
  for (v in c("pre_drybio_ag", "proj_drybio", "T2_drybio_ag")) set(d, which(is.na(d[[v]])), v, 0)
  grm <- rbindlist(lapply(strsplit(o$donors, ",")[[1]], function(s) readRDS(file.path(o$work, sprintf("grm_cond_%s.rds", s)))), fill = TRUE)
  grm <- unique(grm[, .(T2_PLT_CN = as.character(PLT_CN), d_CONDID = as.integer(CONDID), g_tot = gross_tot, m_tot = mort_tot)],
                by = c("T2_PLT_CN", "d_CONDID"))
  d[, d_CONDID := as.integer(d_CONDID)]
  d <- grm[d, on = c("T2_PLT_CN", "d_CONDID")]
  d[, hv := fcoalesce(as.logical(was_harvested), FALSE) & fcoalesce(harvest_intensity, 0) > 0]
  d[, um := fcoalesce(as.logical(was_unmatched), FALSE)]
  d[, h := fifelse(hv, pmin(harvest_intensity, 0.999), 0)]
  d[, net := fifelse(um, 0, proj_drybio / (1 - h) - pre_drybio_ag)]
  d[, base := fifelse(hv, T2_drybio_ag, pre_drybio_ag)]
  d[, r := T2_drybio_ag / d_drybio_ag]
  d[, s_eff := fifelse(is.finite(r) & is.finite(gr_drybio) & abs(r - 1) > 1e-6, (gr_drybio - 1) / (r - 1), 1)]
  d[, donor_ok := !um & !is.na(g_tot) & is.finite(d_drybio_ag) & d_drybio_ag > 0]
  d[, G := fifelse(donor_ok, base * g_tot / d_drybio_ag * s_eff, 0)]
  d[, M := fifelse(donor_ok, base * m_tot / d_drybio_ag * s_eff, 0)]
  d[, resid := net - (G - M)]
  # cem_v3 F3 (track E, docs/PERSEUS_ME_V23_E_PREREG.md amendment 1): on rows with a donor GRM row and no harvest,
  # mortality = donor GRM mortality + the engine's reductions after the donor ratio (eng < 0: SDImax cap, productivity
  # multiplier < 1, disturbance module); gross = net + mortality, so the donor stock-versus-GRM gap goes to gross with
  # its sign. Other rows keep the s41d8 rule.
  d[, ref := fifelse(!um & is.finite(gr_drybio), pre_drybio_ag * (gr_drybio - 1), 0)]
  d[, eng := net - ref]
  d[, v3 := donor_ok & !hv]
  d[, mort := fifelse(v3, M + pmax(-eng, 0), M + pmax(-resid, 0))]
  d[, gross := fifelse(v3, net + mort, G + pmax(resid, 0))]
  d[, rem := fifelse(hv, proj_drybio / (1 - h) * h, 0)]
  d[, loc_key := paste(as.integer(STATECD), as.integer(COUNTYCD), as.integer(PLOT), sep = "_")]
  d[, CONDID := as.integer(CONDID)]
  # per condition x sim x cycle annual rates (cycle = 5 yr); stock at start of cycle 1
  cs <- d[, .(stock = pre_drybio_ag[1], net = sum(net) / 5, gross = sum(gross) / 5, mort = sum(mort) / 5,
              rem = sum(rem) / 5, end = sum(proj_drybio)), by = .(loc_key, CONDID, sim, cycle)]
  cs
}
cs <- tool(o$archive)
sims <- sort(unique(cs$sim)); nS <- length(sims)

# FIA expansion: EXPVOL evaluation of the eval year, condition adjustment mix as in the scorer
plt <- fread(f(ST, "PLOT"), select = c("CN", "STATECD", "COUNTYCD", "PLOT"), colClasses = list(character = "CN"))
plt[, loc_key := paste(STATECD, COUNTYCD, PLOT, sep = "_")]
cnd <- fread(f(ST, "COND"), select = c("PLT_CN", "CONDID", "COND_STATUS_CD"), colClasses = list(character = "PLT_CN"))
peg <- fread(f(ST, "POP_EVAL_GRP"), colClasses = list(character = "CN"))
pe <- fread(f(ST, "POP_EVAL"), colClasses = list(character = c("CN", "EVAL_GRP_CN")))
pet <- fread(f(ST, "POP_EVAL_TYP"), colClasses = list(character = c("EVAL_CN", "EVAL_GRP_CN")))
ev <- merge(pet[, .(EVAL_CN, EVAL_TYP)], pe[EVAL_GRP_CN %in% peg[EVAL_GRP == sc * 10000 + as.integer(o$eval_year), CN], .(CN, EVALID)],
            by.x = "EVAL_CN", by.y = "CN")
ps <- fread(f(ST, "POP_STRATUM"), colClasses = list(character = "CN"))
ppsa <- fread(f(ST, "POP_PLOT_STRATUM_ASSGN"), select = c("PLT_CN", "STRATUM_CN", "EVALID"), colClasses = list(character = c("PLT_CN", "STRATUM_CN")))
id <- ev[EVAL_TYP == "EXPVOL", EVALID][1]
wv <- unique(ppsa[EVALID == id][ps[, .(STRATUM_CN = CN, EXPNS, A1 = ADJ_FACTOR_SUBP, A2 = ADJ_FACTOR_MICR, A3 = ADJ_FACTOR_MACR)],
                                  on = "STRATUM_CN", nomatch = NULL], by = "PLT_CN")
wv <- plt[, .(PLT_CN = CN, loc_key, COUNTYCD)][wv, on = "PLT_CN"]
if (length(cty)) wv <- wv[COUNTYCD %in% cty]
stk <- readRDS(file.path(o$work, sprintf("stock_cond_%s.rds", ST)))
ec <- merge(cnd[COND_STATUS_CD == 1, .(PLT_CN, CONDID = as.integer(CONDID))], wv, by = "PLT_CN")
ec <- merge(ec, stk[, .(PLT_CN, CONDID = as.integer(CONDID), s1, s2, s3)], by = c("PLT_CN", "CONDID"), all.x = TRUE)
for (v in c("s1", "s2", "s3")) set(ec, which(is.na(ec[[v]])), v, 0)
ec[, adj := fifelse(s1 + s2 + s3 > 0, (s1 * A1 + s2 * A2 + s3 * A3) / (s1 + s2 + s3), A1)]
Q <- c("stock", "net", "gross", "mort", "rem", "end")
expand <- function(cs) { nS <- uniqueN(cs$sim)
  tv <- merge(ec[, .(PLT_CN, EXPNS, loc_key, CONDID, adj)], cs, by = c("loc_key", "CONDID"), allow.cartesian = TRUE)
  tv[, lapply(.SD, function(x) sum(x * adj * EXPNS) * LB_TG / nS), by = cycle, .SDcols = Q][order(cycle)] }
tot <- expand(cs)
cov <- c(frame_plots = uniqueN(wv$PLT_CN), forest_conditions = nrow(ec),
         area_share_in_tool = ec[, sum(EXPNS * (paste(loc_key, CONDID) %in% paste(cs$loc_key, cs$CONDID))) / sum(EXPNS)])
print(tot); print(cov)

# contract rows: start = cycle 1 stock; period rows per 5 yr cycle; gate = mean of cycles 1 and 2 (scorer: sum / 10)
g12 <- tot[cycle %in% 1:2]
rows <- rbind(
  data.table(row_type = "start", year = 0L, period_years = 0L, stock = tot[cycle == 1, stock]),
  tot[, .(row_type = "period", year = 5L * cycle, period_years = 5L, stock = end, gross_growth = gross, mortality = mort,
          removals = rem, net_growth = net)],
  data.table(row_type = "gate", year = 10L, period_years = 10L, stock = tot[cycle == 1, stock], gross_growth = mean(g12$gross),
             mortality = mean(g12$mort), removals = mean(g12$rem), net_growth = mean(g12$net),
             note = sprintf("registered s41d8 version A; BAU arm, removals added back; frame plots %d, area share in tool %.4f",
                            cov[["frame_plots"]], cov[["area_share_in_tool"]])), fill = TRUE)
rows[, `:=`(mortality_definition = "cem_v3: donor GRM mortality scaled to subject plus engine reductions after the donor ratio; s41d8 on harvested rows",
            unit = "Tg dry biomass; rates per yr", foliage_in_stock = FALSE, stock_excl_foliage = stock, standing_dead = NA_real_)]
rows[row_type == "period", stock := NA_real_]  # end of cycle live stock is not net of harvest; reported in trajectory
P <- c(engine = "cem", engine_version = o$engine_version, arm = o$arm, aoi_id = o$aoi_id, config_md5 = o$config_md5, run_id = o$run_id)
G <- c(names(P), "row_type", "year", "period_years", "stock", "stock_excl_foliage", "gross_growth", "mortality",
       "mortality_definition", "removals", "net_growth", "standing_dead", "unit", "foliage_in_stock", "note")
for (k in names(P)) rows[[k]] <- P[[k]]
for (k in setdiff(G, names(rows))) rows[[k]] <- NA
dir.create(o$out, showWarnings = FALSE, recursive = TRUE)
fwrite(rows[, ..G], file.path(o$out, "gate_components.csv"))
trj <- function(tt, scen, clim) rbind(data.table(year = 0L, value = tt[cycle == 1, stock]), tt[, .(year = 5L * cycle, value = end)])[, `:=`(scenario = scen, climate = clim)]
tr <- trj(tot, "gate_bau", "rcp45")
if (!is.null(o$scen_archives)) for (z in strsplit(o$scen_archives, ";")[[1]]) {  # scenario:climate=path
  k <- strsplit(z, "=")[[1]]; sc2 <- strsplit(k[1], ":")[[1]]; tr <- rbind(tr, trj(expand(tool(k[2])), sc2[1], sc2[2])) }
tr <- tr[, .(engine = "cem", engine_version = o$engine_version, arm = o$arm, aoi_id = o$aoi_id, config_md5 = o$config_md5,
             run_id = o$run_id, scenario, climate, year, pool = "ag_live", value, unit = "Tg dry biomass")]
fwrite(tr, file.path(o$out, "trajectory.csv"))
jsonlite::write_json(c(as.list(P), list(archive = o$archive, archive_md5 = unname(tools::md5sum(o$archive)), counties = o$counties,
  eval_year = o$eval_year, donors = o$donors, n_sims = nS, coverage = as.list(cov))), file.path(o$out, "engine_meta.json"),
  auto_unbox = TRUE, pretty = TRUE)
cat("CONVERT_CEM_DONE\n")
