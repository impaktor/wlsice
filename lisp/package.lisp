;;; Package definition

(in-package :common-lisp-user)

(defpackage :weighted-least-squares-including-correlation-in-error
  (:use :common-lisp)
  (:nicknames :wls-ice)
  (:export fit))
