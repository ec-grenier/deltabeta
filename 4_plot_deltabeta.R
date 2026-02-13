#==============================================================================#
#' Plot Mortality Curves
#'
#' to run: 
#'   - set regions list
#'   - set list of scenarios
#==============================================================================#

#==============================================================================#
# 0. Source, set packages and paths -----------

packages = c("data.table", "dplyr", "tidyr", "ggplot2", "glue", "cowplot", "purrr", "gridExtra")
invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))
rm(packages)

# I/O
REPO = '/project/cil/home_dirs/egrenier/cil-comms/deltabeta'
input = glue('{REPO}/output/db_out')
out = glue('{REPO}/output/figures')

#==============================================================================#
# 1. Set run specs -----------

# ---- enter region list ---- # 
region_list = c('CAN.11.259.4274', 'CAN.11.269.4448')
#region_list = fread(glue('{REPO}/data/misc/hierarchy-flat.csv'))$`region-key`[2200:2205] # can pick anything in [1, 24378]

# ---- set specs---- # 
# available scenarios:
#  - agegroup: young, older, oldest, combined
#  - rcp:      rcp45, rcp85
#  - gcm:      CCSM4
#  - ssp:      SSP2, SSP3
#  - year:     2025, 2050, 2099
specs = expand.grid(agegroup = c('young', 'older', 'oldest', 'combined'), 
                    rcp = c('rcp85'), 
                    gcm = c('CCSM4'),
                    iam = c('low'),
                    ssp = c( 'SSP3'),
                    year = c(2099)) 

#==============================================================================#
# 2. Plotting functions ----

plot_curve = function(df, year,
                      x.lim=c(-20,45),
                      y.lim=c(-5,200),
                      colors = rev(c("#ff6961","#FBC17D", "#81176D")),
                      y.lab="Change in deaths / 100,000",
                      c.margin=c(.25,.25,0,.25)) {
  
  # Extract numeric temperature from bin, get long df
  df[, temp := as.numeric(sub("\\((-?\\d+),.*", "\\1", bin)) + 0.5]
  df_long = melt(df, id.vars = c("region", "bin", "temp"),
                 measure.vars = c("beta_fa", "beta_ia", "beta_na"),
                 variable.name = "beta", value.name = "value")
  
  # Map beta names to year labels
  df_long[, year := fcase(
    beta == "beta_fa", "Full Adaptation",
    beta == "beta_ia", "Income Adaptation",
    beta == "beta_na", "No Adaptation"
  )]
  
  p = ggplot(data = df_long, aes(x = temp, y = value, group = year)) +
    geom_line(aes(colour = factor(year)), linewidth = 1) +
    #geom_step(aes(colour = factor(year)), linewidth = 1) + # for step-curves
    geom_hline(yintercept = 0, linewidth = .2) +
    theme_minimal() +
    ylab(y.lab) +
    coord_cartesian(ylim = y.lim, xlim = x.lim) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = y.lim) +
    scale_linetype_discrete(name = NULL) +
    scale_color_manual(values = colors, name = "Curve type") +
    theme(legend.justification = c(0, 1),
          legend.position = c(0.05, 0.95),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          plot.margin = unit(c.margin, "in"))
  
  return(p)
  
}

plot_hist = function(df, year,
                     x.lim=c(-20,45), 
                     hist.y.lim=c(0,50), 
                     hist.breaks = seq(-20, 50, by = 10), 
                     h.margin=c(.05,.25,.25,.25),
                     hist.x.lab = "Daily temperature (C)",
                     hist.y.lab = "Number of days",
                     colors = rev(c("#ff6961","#FBC17D", "#81176D"))) {

    t_col = paste0("T_", year)
    
    # Filter out summary rows and melt to long format
    df = df[!bin %in% c("Total <20C", "Total >20C", "Total")]
    df_long = melt(df, id.vars = c("region", "bin"),
                   measure.vars = c("T_1993", "T_2005", t_col),
                   variable.name = "year", value.name = "value")
    
    # Clean up year labels, get midpoint for each bar
    df_long[, year := gsub("T_", "", year)]
    df_long[, temp := as.numeric(sub("\\((-?\\d+),.*", "\\1", bin)) + 0.5]
    
    # plot
    p = ggplot(data = df_long) +
      geom_bar(aes(x = temp, y = value, fill = factor(year)), 
               stat = "identity", position = "dodge", alpha = 1, orientation = "x") +
      theme_minimal() +
      scale_x_continuous(expand = c(0, 0)) + 
      scale_y_continuous(expand = c(0, 0), breaks = hist.breaks) + 
      ylab(hist.y.lab) + xlab(hist.x.lab) +
      coord_cartesian(xlim = x.lim, ylim = hist.y.lim) + 
      scale_fill_manual(values = colors, name = "Year") +
      theme(legend.justification = c(0, 1),
            legend.position = c(0.05, 0.95),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            plot.margin = unit(h.margin, "in"))
    
    return(p)
}

plot_db = function(data, reg, agegroup, year, rcp, gcm, iam, ssp, out) {

  message(glue('[running: ] {reg}, {year}, {rcp}, {gcm}, {iam}, {ssp}'))
  
  # filter inputted data
  data = data[data$region == reg, ]
  db = data[!bin %in% c("Total <20C", "Total >20C", "Total")]
  summary = data[bin %in% c("Total <20C", "Total >20C", "Total"), .(bin, effect_fa, effect_ia, effect_na)]
  setnames(summary, 
           c("bin", "effect_fa", "effect_ia", "effect_na"),
           c("", "Full Adaptation", "Income Adaptation", "No Adaptation"))
  
  # create curve and histogram plot
  message('\nplotting...\n')
  
  curve = plot_curve(df=db, year=year)
  hist = plot_hist(df=db, year=year)
  summary = tableGrob(summary, rows = NULL, theme = ttheme_default(base_size = 10))
  
  # generate title, combine plots
  title = ggdraw() + draw_label(glue('{reg} {year} {rcp} {gcm} {iam} {ssp} {agegroup}'), fontface = 'bold')
  deltabeta = plot_grid(title, curve, hist, summary, ncol = 1, rel_heights = c(0.1, 1, 0.5, 0.4))
  
  message('done.\n')
  print(deltabeta)
  
  # save
  message(glue('\n==== saving ==== \n\nOutput directory: {out} \n'))
  message(glue('File name: {reg}-delta_beta-{year}-{rcp}-{gcm}-{iam}-{ssp}-{agegroup}-plot.pdf'))

  fname = glue('{reg}-delta_beta-{year}-{rcp}-{gcm}-{iam}-{ssp}-{agegroup}-plot.pdf')
  ggsave(filename=fname, path=out, plot=deltabeta, height=10, width=8)
  
  message('\n==== saved ====\n')
  
}

#==============================================================================#
# 3. Run plots ----

# If running plots for multiple regions and one scenario, data will get loaded in only once 
pwalk(specs, function(agegroup, year, rcp, gcm, iam, ssp) {
  
  data = fread(glue('{input}/{rcp}/{gcm}/{iam}/{ssp}/mortality-delta_beta-{year}-{rcp}-{gcm}-{iam}-{ssp}-{agegroup}.csv'))
  
  for (region in region_list) {
    plot_db(data, region, agegroup, year, rcp, gcm, iam, ssp, out)
  }
  
})

