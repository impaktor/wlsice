;;;; Least squares including correlation in error

(in-package ls-ice)

(defun compute-mean (trajectories)
  "Return foreign array of dim. 1xN of MxN input,
 mean over the M samples for each N."
  (let ((mean-of-trajectories '())
        (N-sampling-points (cadr (grid:dimensions trajectories))))
    (dotimes (i N-sampling-points (nreverse mean-of-trajectories))
      (push (gsll:mean (gsll:column trajectories i)) mean-of-trajectories))))

(defvar matrix2 (GRID:MAKE-FOREIGN-ARRAY 'DOUBLE-FLOAT :INITIAL-CONTENTS
                                         '((-43.5d0 4.2d0 1.2d0)
                                           (-83d0 3.12d0 -6.1d0)
                                           (94.72d0 31.94d0 3.2d0))))


;; TO-DO: test it!
(defun make-covariance-matrix (trajectories)
  "Compute the covariance matrix estimator, from a foreign array,
 (in a very inefficient way)"
  (let ((N (cadr (grid:dimensions trajectories))))
    (loop for i from 0 below N
       collect (loop for j from 0 below N
                  collect (gsll:covariance (gsll:column trajectories i)
                                           (gsll:column trajectories j))))))


(defun compute-ls-ice-error (df d2f R C)
  "Compute the correlation corrected error estimate")


(defun fit-function (time y-mean R C)
  "Fit the data to the function")


(defun fit (f df d2f time trajectories )
  "the exported function to be used"  )


(defun test () (print "test 42!"))

(defvar *data* "../data/" "Path to the data files")
