;;;; wls-ice.asd

;;; small nice example: https://github.com/BradWBeer/indy-anna

(asdf:defsystem #:weigted-least-squares-including-correlation-in-error
  :description "Applies the weighted least squares, including correlation in error, to data of ensamble averages for correct error estimation"
  :author "Karl Fogelmark"
  :license "GPL"
  :depends-on (#:gsll)
  :serial t
  :components ((:file "wls-ice")
               (:file "example")))
