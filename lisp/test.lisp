;;; My file with test code.

(asdf:operate 'asdf:load-op :gsll)

;; How to look up documentation from GSL doc to GSLL
(gsl:gsl-lookup "gsl_matrix_get_col") ;; -> GSLL:COLUMN
(documentation #'GSLL:COLUMN 'function)

;; Here is a 3x3 GSLL array:
(defvar matrix1 (GRID:MAKE-FOREIGN-ARRAY 'DOUBLE-FLOAT :INITIAL-CONTENTS
                                        '((-34.5d0 8.24d0 3.29d0)
                                          (-8.93d0 34.12d0 -6.15d0)
                                          (49.27d0 -13.49d0 32.5d0))))

(defvar matrix2 (GRID:MAKE-FOREIGN-ARRAY 'DOUBLE-FLOAT :INITIAL-CONTENTS
                                        '((-43.5d0 4.2d0 1.2d0)
                                          (-83d0 3.12d0 -6.1d0)
                                          (94.72d0 31.94d0 3.2d0))))

(grid:dimensions matrix1)
(grid:gref matrix1 0 0)


;;========================================

; (mapcar (lambda (x) (format t "~s~%" x)) (gsl:examples))


;; Why the "copy to"?
(GRID:COPY-TO (GSLL:MATRIX-TRANSPOSE GSLL::M1 GSLL::M2))


;; (in-package antik-user)

;; (setf grid:*default-grid-type* 'grid:foreign-array)
;; (grid:sin (GRID:MAKE-FOREIGN-ARRAY 'DOUBLE-FLOAT :INITIAL-CONTENTS
;;                               '((-34.5d0 8.24d0 3.29d0)
;;                                 (-8.93d0 34.12d0 -6.15d0)
;;                                 (49.27d0 -13.49d0 32.5d0))))

;; (invert-matrix
;;    (grid:make-foreign-array
;;     'double-float :initial-contents '((1.0d0 2.0d0) (3.0d0 4.0d0))))

;; Om(?) Antik:
;; file:///home/vandelay/quicklisp/dists/quicklisp/software/antik-master-ad6432e3-git/documentation/build/html/gridintro.html
;; The argument grid-type to a function or the variable
;; grid:*default-grid-type* should be bound to one of 'cl:array or
;; 'grid:foreign-array, to specify an array or foreign-array
;; respectively.



;; ;;; makes matrix m1 zero
;; (let ((gsll::m1
;;         (grid:make-foreign-array 'double-float :initial-contents
;;                                  '((-34.5d0 8.24d0 3.29d0)
;;                                    (-8.93d0 34.12d0 -6.15d0)
;;                                    (49.27d0 -13.49d0 32.5d0)))))
;;    (gsll:set-zero gsll::m1)
;;    (grid:copy-to gsll::m1))


;; COV
(LET ((GSLL::V1
       (GRID:MAKE-FOREIGN-ARRAY 'DOUBLE-FLOAT :INITIAL-CONTENTS
                                '(-34.5d0 8.24d0 3.29d0)))
      (GSLL::V2
       (GRID:MAKE-FOREIGN-ARRAY 'DOUBLE-FLOAT :INITIAL-CONTENTS
                                '(-8.93d0 34.12d0 -6.15d0))))
  (LET ((GSLL::MEAN1 (GSLL:MEAN GSLL::V1))
        (GSLL::MEAN2 (GSLL:MEAN GSLL::V2)))
    (GSLL:COVARIANCE GSLL::V1 GSLL::V2 GSLL::MEAN1 GSLL::MEAN2)))

;; GSLL::NONLINEAR-LEAST-SQUARES
;; GSLL::LINEAR-LEAST-SQUARES
;; GSLL::MINIMIZATION-MULTI
;; GSLL::MINIMIZATION-ONE

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
