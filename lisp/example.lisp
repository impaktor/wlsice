;;; Example file, that reads in data from folder ../data, and fits a non-linear function: a*t**b to it.

(in-package :common-lisp-user)

(defun f (t parameters)
  "The analytical function of t to fit parameters to"
  (let ((a (car parameters))
        (b (cadr parameters)))
    (mapcar (lambda (x) (* a (expt x b))) t)))

(defun df (t parameters)
  "Gradient of analytical function of t. Return 2-element list of gradient in each point t for each of our two parameters a, b"
  (let ((a (car parameters))
        (b (cadr parameters)))
    (list (mapcar (lambda (x) (expt x b)) t)
          (mapcar (lambda (x) (* a (log x) (expt x b)))))))

(defun d2f (t parameters)
  "Hessian matrix of analytical function to fit"
  (let ((a (car parameters))
        (b (cadr parameters)))
    ;; TODO!
    (mapcar (lambda (x) 0) t)           ; 0,0 is all 0
    (mapcar (lambda (x) (* (log x) (expt x b))) t)    ;0,1
    (mapcar (lambda (x) (* (log x) (expt x b))) t)    ;1,0
    (mapcar (lambda (x) (* a (expt (log x) 2) (expt x b))))  ;1,1
    ))

(defun read-in-data (path)
  "Reads in two-column data from path, returns a suitable data structure, of time vector and matrix of trajectories")

;; start-time 200

;(ccls:fit f df d2f time trajectories :guess '(1 1))
