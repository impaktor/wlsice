; (asdf:operate 'asdf:load-op :gsll)
; (mapcar (lambda (x) (format t "~s~%" x)) (gsl:examples))

;;; makes matrix m1 zero
(let ((gsll::m1
        (grid:make-foreign-array 'double-float :initial-contents
                                 '((-34.5d0 8.24d0 3.29d0)
                                   (-8.93d0 34.12d0 -6.15d0)
                                   (49.27d0 -13.49d0 32.5d0)))))
   (gsll:set-zero gsll::m1)
   (grid:copy-to gsll::m1))

;; GSLL::NONLINEAR-LEAST-SQUARES
;; GSLL::LINEAR-LEAST-SQUARES
;; GSLL::MINIMIZATION-MULTI
;; GSLL::MINIMIZATION-ONE
;; GSLL:COVARIANCE
;; GSLL:AUTOCORRELATION
;; GSLL::CHI-SQUARED
;; GSLL::EXPONENTIAL-POWER
;; GSLL::EXPONENTIAL
;; GSLL:CDOT
;; GSLL::DOT
;; GSLL::POWER
;; GSLL::LOGARITHM
;; GSLL::EXPONENTIAL-FUNCTIONS
;; GSLL::POLYNOMIAL
;; GSLL:MATRIX-TRANSPOSE
;; GSLL:MATRIX-TRANSPOSE*
;; GSLL::SETF-COLUMN
;; GSLL:COLUMN
;; GSLL::SETF-ROW
;; GSLL:ROW
;; GSLL:SET-IDENTITY
;; GSLL::VECTOR-SET-ZERO
;; GSLL::MATRIX-SET-ALL
;; GSLL::VECTOR-SET-ALL
