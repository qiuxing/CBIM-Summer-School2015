(TeX-add-style-hook "article"
 (lambda ()
    (LaTeX-add-environments
     "kframe"
     "knitrout")
    (TeX-add-symbols
     '("hlkwd" 1)
     '("hlkwc" 1)
     '("hlkwb" 1)
     '("hlkwa" 1)
     '("hlstd" 1)
     '("hlopt" 1)
     '("hlcom" 1)
     '("hlstr" 1)
     '("hlnum" 1)
     "maxwidth"
     "FrameCommand")
    (TeX-run-style-hooks
     "alltt"
     "framed"
     "color"
     "hyperref"
     "xspace"
     "calc"
     "ifthen"
     "bm"
     "subfigure"
     "graphicx"
     "mathrsfs"
     "amsthm"
     "amssymb"
     "amsmath"
     "amsfonts"
     "wordlike"
     "fontenc"
     "T1"
     "times"
     "inputenc"
     "latin1"
     "babel"
     "english"
     "beamerarticle"
     "latex2e"
     "art10"
     ""
     "main")))

