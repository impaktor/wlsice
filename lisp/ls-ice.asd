;;;; ls-ice.asd

;;; small nice example: https://github.com/BradWBeer/indy-anna

(asdf:defsystem #:least-squares-including-correlation-in-error
  :description "Applies the least squares including correlationin error to data of ensamble averages for correct error estimation"
  :author "Karl Fogelmark"
  :license "GPL"
  :depends-on (#:gsll)
  :serial t
  :components ((:file "ls-ice")
               (:file "example")))
