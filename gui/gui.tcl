#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" ${1+"$@"}

package require Tk

# -----------------------------------------------------------------------------
# @brief  GUI for the Schwarzschild case
#
#         Just write the parameters in a parameters.txt file
# -----------------------------------------------------------------------------

# Common data for the calculations
set uiCommon(RangoInf) -60.0
set uiCommon(RangoSup) 60.0
set uiCommon(l) 0.0
set uiCommon(rS) 2.0
set uiCommon(wMin) 0.02
set uiCommon(wMax) 1.0
set uiCommon(nW) 100
set uiCommon(nL) 5
# Common data for the wave graph
set uiCommon(X_LOW) -60.0
set uiCommon(X_HIG) 60.0
set uiCommon(Y_LOW) -2.0
set uiCommon(Y_HIG) 2.0

# -----------------------------------------------------------------------------
# @brief  Write the gnuplot configuration file
# -----------------------------------------------------------------------------
proc writeGnuplotFile {fileName fileContent} {
  set fd [open $fileName w+]
  
  foreach line $fileContent {
    puts $fd $line
  }
  
  close $fd
}

# -----------------------------------------------------------------------------
# @brief  Modify some gnuplot parameters rewritting the gnuplot files
# -----------------------------------------------------------------------------
proc writeGnuplotSettings {} {
  # Tortoise
  set fileContent {}
  lappend fileContent {set title 'Tortoise coordinates'}
  lappend fileContent "set xrange \[$::uiCommon(X_LOW) : $::uiCommon(X_HIG)\]"
  lappend fileContent {set yrange [-5 : 60.0]}
  lappend fileContent {plot 'turtle.txt' using 1:2 with lines lc "blue" title "x[y] calculated", 'turtle.txt' using 3:1 with lines lc "green" title "y[x] ana.inv."}

  writeGnuplotFile "plot-turtle.gnu" $fileContent

  # Wave real / imaginary
  lappend fileContent {set title 'Schwarzschild'}
  lappend fileContent "set xrange \[$::uiCommon(X_LOW) : $::uiCommon(X_HIG)\]"
  lappend fileContent "set yrange \[$::uiCommon(Y_LOW) : $::uiCommon(Y_HIG)\]"
  lappend fileContent {plot 'wave.txt' using 1:2 with lines lc "green" title "Φ Re", 'wave.txt' using 1:3 with lines lc "red" title "Φ Im", 'wave.txt' using 1:4 with lines title "V(x)"}
  
  writeGnuplotFile "plot.gnu" $fileContent

  # Potential
  set fileContent {}
  lappend fileContent {set title 'Potential'}
  lappend fileContent "set xrange \[$::uiCommon(X_LOW) : $::uiCommon(X_HIG)\]"
  lappend fileContent {set yrange [-0.05 : 1.2]}
  lappend fileContent {plot 'wave.txt' using 1:4 with lines lc "blue" title "V(x)"}
  
  writeGnuplotFile "plot-v.gnu" $fileContent
  
  # RT coeffcients
  set fileContent {}
  lappend fileContent {set title 'Coeficientes R y T'}
  lappend fileContent {set xrange [0.0 : 1.0]}
  lappend fileContent {set yrange [0.0 : 1.2]}
  lappend fileContent {set xlabel 'ω'}
  lappend fileContent {plot 'coefficients.txt' using 1:2 with lines lc "green" title "R", 'coefficients.txt' using 1:3 with lines lc "red" title "T", 'coefficients.txt' using 1:4 with lines title "R + T"}
  
  writeGnuplotFile "plot-k.gnu" $fileContent
  
  # Error's Log10 for R + T = 1
  set fileContent {}
  lappend fileContent {set title 'log10|1 - R - T|'}
  lappend fileContent "set xrange \[$::uiCommon(wMin) : $::uiCommon(wMax)\]"
  lappend fileContent {set yrange [-20.0 : 0]}
  lappend fileContent {set xlabel 'ω'}
  lappend fileContent {plot 'coefficients.txt' using 1:5 with lines lc "red" title "log10 Err"}
  
  writeGnuplotFile "plot-k-err.gnu" $fileContent
  
  # Sigma
  set fileContent {}
  lappend fileContent {set title 'σl'}
  lappend fileContent {set xrange [0.0 : 1.0]}
  lappend fileContent {set yrange [0.0 : 120.0]}
  lappend fileContent {set xlabel 'ω'}
  lappend fileContent {plot 'sigma-l.txt' using 1:2 with lines lc "red" title "σl"}
  
  writeGnuplotFile "plot-sigma.gnu" $fileContent
}

# -----------------------------------------------------------------------------
# @brief  Executables with options
# -----------------------------------------------------------------------------
proc runExe {fileNameAndOption {delay 3000}} {
  after $delay
  exec {*}$fileNameAndOption &
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Commands
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# -----------------------------------------------------------------------------
# @brief  Save the parameters in the 'parameters.txt' file, do the calculations
#         using the exe and execute gnuplot to visualize the data.
#
#         The order and number of parameters have to be the same as in 'ui.c' file
# -----------------------------------------------------------------------------
proc cmd_visualizar {} {
  set fd [open parameters.txt w+]
  
  puts $fd $::uiCommon(RangoInf); #0
  puts $fd $::uiCommon(RangoSup); #1
  puts $fd $::uiCommon(l);        #2
  puts $fd $::uiCommon(rS);       #3
  puts $fd $::uiCommon(wMin);     #4
  puts $fd $::uiCommon(wMax);     #5
  puts $fd $::uiCommon(nW);       #6
  puts $fd $::uiCommon(nL);       #7
  
  close $fd
  
  # Do the calculations
  set command [list schwarzschild.exe]
  exec {*}$command
  after 1000

  writeGnuplotSettings
  cmd_replot
}

# -----------------------------------------------------------------------------
# @brief  Redraw the graphics
# -----------------------------------------------------------------------------
proc cmd_replot {} {
  # Tortoise coordinates
  set gnuplotCmd [list gnuplot.exe -p plot-turtle.gnu]
  runExe $gnuplotCmd

  # Potential and wave function
  set gnuplotCmd [list gnuplot.exe -p plot.gnu]
  runExe $gnuplotCmd
  
  set gnuplotCmd [list gnuplot.exe -p plot-v.gnu]
  runExe $gnuplotCmd

  # RT coefficients
  set gnuplotCmd [list gnuplot.exe -p plot-k.gnu]
  runExe $gnuplotCmd

  # Log10 err RT
  set gnuplotCmd [list gnuplot.exe -p plot-k-err.gnu]
  runExe $gnuplotCmd

  # Sigma
  set gnuplotCmd [list gnuplot.exe -p plot-sigma.gnu]
  runExe $gnuplotCmd
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# UI
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# -----------------------------------------------------------------------------
# @brief  Modify the numerical methods parameters
# -----------------------------------------------------------------------------
proc ui_numerico {} {
  set mainFrame [frame .ntb.frmNumerico.frmMain]

  set w1 [label $mainFrame.lblRango -text "Range: "]
  set w2 [label $mainFrame.lblRangoInf -text "Inf"]
  set w3 [entry $mainFrame.entRangoInf -textvariable uiCommon(RangoInf) -width 5 -background white]
  set w4 [label $mainFrame.lblRangoSup -text "Sup"]
  set w5 [entry $mainFrame.entRangoSup -textvariable uiCommon(RangoSup) -width 5 -background white]
  grid $w1 $w2 $w3 $w4 $w5
    
  set w1 [button $mainFrame.btnVisualizar -text "Show" -command cmd_visualizar]
  set w2 [button $mainFrame.btnRePlot -text "Re Plot" -command cmd_replot]
  grid $w1 $w2
  
  # Layout
  grid $mainFrame -sticky ew
}

# -----------------------------------------------------------------------------
# @brief  Modify some aspects of the gnuplot graphic
# -----------------------------------------------------------------------------
proc ui_fisico {} {
  set mainFrame [frame .ntb.frmFisico.frmMain]
  
  foreach id {l nL rS wMin wMax nW} {
    set w1 [label $mainFrame.lbl$id -text "$id: "]
    set w2 [entry $mainFrame.ent$id -textvariable uiCommon($id) -width 5 -background white]
    
    grid $w1 $w2
  }
  
  # Layout
  grid $mainFrame -sticky ew
}

# -----------------------------------------------------------------------------
# @brief  Modify some phisical parameters
# -----------------------------------------------------------------------------
proc ui_gnuplot {} {
  set mainFrame [frame .ntb.frmGnuplot.frmMain]
  
  foreach id {X_LOW X_HIG Y_LOW Y_HIG} {
    set w1 [label $mainFrame.lbl$id -text "$id: "]
    set w2 [entry $mainFrame.ent$id -textvariable uiCommon($id) -width 5 -background white]
    
    grid $w1 $w2
  }
  
  # Layout
  grid $mainFrame -sticky ew
}

# -----------------------------------------------------------------------------
# @brief  Init the GUI
# -----------------------------------------------------------------------------
proc ui_init {} {
  ttk::notebook .ntb
  .ntb add [frame .ntb.frmNumerico] -text "Numerical"
  .ntb add [frame .ntb.frmFisico] -text "Physic"
  .ntb add [frame .ntb.frmGnuplot] -text "gnuplot"
  .ntb select .ntb.frmNumerico
  
  ui_numerico
  ui_fisico
  ui_gnuplot

  grid .ntb
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Main
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc main {} {
  ui_init
}

main
