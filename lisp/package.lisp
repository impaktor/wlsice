;;; Package definition

(in-package :common-lisp-user)

(defpackage :least-squares-including-correlation-in-error
  (:use :common-lisp)
  (:nicknames :ls-ice)
  (:export fit))
