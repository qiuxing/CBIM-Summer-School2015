(TeX-add-style-hook "main"
 (lambda ()
    (LaTeX-add-bibliographies
     "xing")
    (LaTeX-add-labels
     "sec:linear-models"
     "eq:lin-mod"
     "eq:Ftest"
     "sec:logistic-models"
     "eq:logistic"
     "sec:categ-data-analys"
     "sec:other-topics")
    (TeX-add-symbols
     '("uvec" ["argument"] 1)
     '("umark" 1)
     '("urepeat" 2)
     '("uinner" 2)
     '("ustab" 1)
     '("uo" 2)
     '("us" 1)
     '("ucorr" 2)
     '("ucor" 2)
     '("ucov" 2)
     '("SL" 1)
     '("GL" 1)
     '("uder" 2)
     '("ensuretext" 1)
     "C"
     "CP"
     "Q"
     "R"
     "RP"
     "Z"
     "med"
     "uPr"
     "uE"
     "uvar"
     "ud"
     "uProj"
     "uimply"
     "uequiv"
     "uforall"
     "ie"
     "eg")))

