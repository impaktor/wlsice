;;;; ccls.asd

;;; small nice example: https://github.com/BradWBeer/indy-anna

(asdf:defsystem #:correlation-corrected-least-squares
  :description "Applies the correlation corrected least squares to data of ensamble averages for correct error estimation"
  :author "Karl Fogelmark"
  :license "GPL"
  :depends-on (#:gsll)
  :serial t
  :components ((:file "ccls")
               (:file "example")))
