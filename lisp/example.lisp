;;; Example file, that reads in data from folder ../data, and fits a non-linear function: a*t**b to it.

(in-package :common-lisp-user)

(asdf:operate 'asdf:load-op :gsll)

;;; Path to trajectory data files
(defparameter *data-path* "../data/")

(defun f (time parameters)
  "The analytical function of t to fit parameters to"
  (let ((a (car parameters))
        (b (cadr parameters)))
    (mapcar (lambda (x) (* a (expt x b))) time)))

(defun df (time parameters)
  "Gradient of analytical function of t. Return 2-element list of gradient
 in each point t for each of our two parameters a, b"
  (let ((a (car parameters))
        (b (cadr parameters)))
    (list (mapcar (lambda (x) (expt x b)) time)
          (mapcar (lambda (x) (* a (log x) (expt x b))) time))))

(defun d2f (time parameters)
  "Hessian matrix of analytical function to fit"
  (let ((a (car parameters))
        (b (cadr parameters)))
    ;; TODO!
    (mapcar (lambda (x) 0) time)           ; 0,0 is all 0
    (mapcar (lambda (x) (* (log x) (expt x b))) time)    ;0,1
    (mapcar (lambda (x) (* (log x) (expt x b))) time)    ;1,0
    (mapcar (lambda (x) (* a (expt (log x) 2) (expt x b))) time)  ;1,1
    ))


(defun get-second-column (input-file)
    "Extract the second column, as a list, from a two-column file with
floats. Skipp the first line in the input-file, it's a comment."
  (with-open-file (stream input-file)
    (read-line stream)                ; first line is a comment, skip it
    (do ((line (read-line stream nil)
               (read-line stream nil))
         (out '()))
        ((null line) (nreverse out))
      (with-input-from-string (in line)
        (read in)                           ; skip first element/column,
        (setf out (cons (read in) out)))))) ; save second element/column

(defun get-trajectories ()
  "Read in directory, and return an MxN matrix of all M trajectories,
each with N sampling points, as a foreign array."
  (grid:make-foreign-array
   'double-float :initial-contents
   (mapcar #'get-second-column
           (directory (concatenate 'string *data-path* "trajectory*")))))

(defvar *trajectories* (get-trajectories) "Foreign MxN array of M trajectories")


;;(wls-ice:fit #'f #'df #'d2f time trajectories :guess '(1 1))
