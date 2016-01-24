require(rkwarddev)
local({
  
  # Author names and contact information
  about.info <- rk.XML.about(
    name = "Melting curve analysis",
    author = c(
      person(given = "Stefan", family = "Roediger",
             email = "Stefan.Roediger@b-tu.de", 
             role = c("aut","cre"))),
    about = list(desc = "GUI interface to perfrom melting curve analysis",
                 version = "0.0.1-2", url = "")
  )
  
  ## help page
  plugin.summary <- rk.rkh.summary(
    "Analysis of melting data. The plugin is primarily targetet at the analysis of melting curve data from nucleic acid experiments but might be used in other cases too."
  )
  plugin.usage <- rk.rkh.usage(
    "Chose a data set and method to analyse melting curve data."
  )
  
  # Define dependencies
  dependencies.info <- rk.XML.dependencies(dependencies = list(rkward.min = "0.6.1"), 
                                           package = list(c(name = "MBmca", min = "0.0.3-5")))
  # General settings
  
  # Definition of plot labels and appearance
  generic.plot.options <- rk.plotOptions()
  plot.main <- rk.XML.input(label = "Main title", initial = "Melting curve analysis")
  plot.xlab <- rk.XML.input(label = "Abscissa label", initial = "Temperature")
  plot.ylab <- rk.XML.input(label = "Ordinate label", initial = "-d(RFU)/dT")
  
  var.select <- rk.XML.varselector(label = "Select data")
  selected.x <- rk.XML.varslot(label = "Temperature", source = var.select, types = "number", required = TRUE)
  var.data <- rk.XML.varslot(label = "Signal", source = var.select, multi = TRUE, classes = "numeric", types = "number", required = TRUE)
  
  
  # Plot preview
  preview.chk <- rk.XML.preview(label = "Preview")
  generic.plot.options <- rk.plotOptions()
  plot.text.color <- rk.plotOptions(embed = "rkward::color_chooser", button = FALSE)
  
  basic.settings <- rk.XML.row(
    var.select,
    rk.XML.col(
      selected.x,
      var.data,
      preview.chk,
      rk.XML.stretch()
    ))
  
  # Defintion of setting for the analysis
  # Definion for the mcaSmoother function
  minmax <- rk.XML.cbox(label = "Min-Max normalization", val = ", minmax  =  TRUE", chk = TRUE)
  df.fact.spin <- rk.XML.spinbox(label = "df.fact", min = "0.6", max = "1", initial = "0.99")
  overplot.raw <- rk.XML.cbox(label = "Show raw data", value = "1", un.value = "0")
  bgadj <- rk.XML.cbox(label = "Adjust background", value = "TRUE", un.value = "FALSE")
  
  # Definion for the diffQ2 function
  inder.chk <- rk.XML.cbox("Use inder function", value = "1", un.value = "0")
  rsm.chk  <- rk.XML.cbox("Double temperature resolution", value = "1", un.value = "0")
  fws.spin <- rk.XML.spinbox(label = "Windowsize for Tm calculation", min = "2", max = "8", initial = "8", real = FALSE)
  fct.drop <- rk.XML.dropdown(label = "Test for minimum, maximum or both",
                              options = list("Minimum" = c(val = "min"), 
                                             "Maximum" = c(val = "max", chk=TRUE),
					     "Both" = c(val = "both")))
  
  legend.pos.drop <- rk.XML.dropdown(label = "Position of legend",
                                     options = list("Bottomright" = c(val = "bottomright"), 
                                                    "Bottom" = c(val = "bottom"),
                                                    "Bottomleft" = c(val = "bottomleft"),
                                                    "Left" = c(val = "left"),
                                                    "Topleft" = c(val = "topleft", chk = TRUE),
                                                    "Top" = c(val = "top"),
                                                    "Topright" = c(val = "topright"),
                                                    "Right" = c(val = "right"),
                                                    "Center" = c(val = "center")))
  ncol.legend.spin <- rk.XML.spinbox(label = "Number of columns in legend", min = "1", initial = "1", real = FALSE)
  legend.frame <- rk.XML.frame(legend.pos.drop, ncol.legend.spin, rk.XML.stretch(), label="Legend")

  rsm.chk  <- rk.XML.cbox("Use rsm", value = "1", un.value = "0")
  warn.chk  <- rk.XML.cbox("Surpress warnings", value = "1", un.value = "0")
  abline.Tm.chk  <- rk.XML.cbox("Show Tm", value = "1", un.value = "0")
  
  full.dialog <- rk.XML.dialog(
    label = "Melting Curve Analysis",
    rk.XML.tabbook(tabs = list("Basic settings" = list(basic.settings), 
                               "Options Smoothing" = list(df.fact.spin, minmax, overplot.raw, bgadj), 
                               "Options MCA" = list(fws.spin, inder.chk, warn.chk, fct.drop),
                               "Plot options" = list(generic.plot.options, abline.Tm.chk, 
                                                     legend.frame)))
  )
  
  JS.calc <- rk.paste.JS(
    js.var.data <- rk.JS.vars(var.data, join = ", "), # get selected vars
    echo("xy <- as.matrix(data.frame(rk.list(", selected.x,", ", js.var.data,")))\n\n"),
    
    echo("
	  res.smooth <- cbind(xy[, 1], apply(xy[, -1], 2, function(i) {
	  xy.in <- mcaSmoother(xy[, 1], i, bgadj = ", bgadj,", df.fact = ", df.fact.spin,")[, 2]
	  }))
    \n"),
    
    echo("res.min <- lapply(2L:ncol(res.smooth), function(i) {\n"),
    echo("\tout <- diffQ2(cbind(res.smooth[, 1], res.smooth[, i])"),
    ite(id(inder.chk), echo(", inder = TRUE")),
    echo(", fct = min, fws = ", fws.spin,""),
    echo(", verbose = TRUE, warn = FALSE)})\n"),

    echo("res.max <- lapply(2L:ncol(res.smooth), function(i) {\n"),
    echo("\tout <- diffQ2(cbind(res.smooth[, 1], res.smooth[, i])"),
    ite(id(inder.chk), echo(", inder = TRUE")),
    echo(", fct = max, fws = ", fws.spin,""),
    echo(", verbose = TRUE, warn = FALSE)})\n"),
    
    echo("range.T <- c(range(xy[, 1]))\n"),
    echo("range.fluo <- c(range(xy[, -1]))\n"),
    echo("range.TmD1.fluo <- c(range(unlist(lapply(1L:length(res.min), function(y) {try(res.min[[y]][[\"TmD1\"]][[\"xy\"]][, 2])}))))\n")
  )
  
  JS.print <- rk.paste.JS(
    rk.paste.JS.graph(
      
      echo("par(mfrow = c(1,2))\n"),
      echo("plot(NA, NA, xlim = range.T, ylim = range.fluo)\n"),
      ite(overplot.raw, 
          echo("lapply(2L:ncol(xy), function(y) {try(lines(xy[, 1], xy[, y], col = 2))})\n")),
      echo("lapply(2L:ncol(res.smooth), function(y) {try(lines(res.smooth[, 1], res.smooth[, y]))})\n"),
      ite(abline.Tm.chk,
          echo("lapply(1L:length(res.min), function(y) {try(abline(v = res.min[[y]]$TmD1$Tm))})\n")),
      echo("legend(\"", legend.pos.drop,"\", colnames(res.smooth[, -1]), ncol = ", ncol.legend.spin,")\n"),
      
      echo("plot(NA, NA, xlim = range.T, ylim = range.TmD1.fluo)\n"),
      echo("lapply(1L:length(res.min), function(y) { try(lines(res.min[[y]]$TmD1$xy))})\n"),
      echo("lapply(1L:length(res.min), function(y) { try(lines(res.min[[y]]$Tm1D2$xy))})\n"),
      ite(abline.Tm.chk,
          echo("lapply(1L:length(res.min), function(y) {try(abline(v = res.min[[y]]$TmD1$Tm))})\n")),
      echo("par(mfrow = c(1,1))\n")
      
    ),
    ite("full", rk.paste.JS(
      echo("\nrk.print(res.min)\n"),
      level = 3)
    )
  )
  
  rk.plugin.skeleton(
    about = about.info,
    dependencies = dependencies.info,
    xml = list(dialog = full.dialog),
    js = list(require = "MBmca",
              calculate = JS.calc,
              doPrintout = JS.print),
    rkh = list(plugin.summary, plugin.usage),
    pluginmap = list(
      name = "Melting curve analysis",
      hierarchy = list("analysis", "MCA")),
    load = TRUE,
    overwrite = TRUE,
    show = TRUE
  )
  
})
