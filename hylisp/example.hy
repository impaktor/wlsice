#!/usr/bin/env hy

"""Example file for using WLS-ICE method to fit a non-linear /
powerlaw function to mean ensamble data.

"""

(import [numpy :as np] sys os glob wls-ice)

;; np.array(1), np.array(1) -> np.array(1)
(defn f [t params]
  """A powerlaw: a*t**b"""
  (assert (len (= 2 params)))
  (np.multiply (get params 0) (np.power t (get params 1))))


;; np.array(1), np.array(1) -> np.array(2)
(defn df [t params]
  """Gradient of a powerlaw, with regards to parameters a,b:
  gradient(a*t**b) = [t**b, a*log(t)*t**b]
  """
  (assert (> (get t 0) 0))
  (assert (= (len params) 2))

  (setv (, a b) (, (get params 0)
                   (get params 1)))

  (setv dfd_lam (np.zeros (, (len params) (len t)))
        (get dfd_lam 0) (np.power t b)
        (get dfd_lam 1) (np.multiply (* a (np.log t))
                                     (np.power t b)))
  dfd_lam)

;; np.array(1), np.array(1) -> np.array(3!)
(defn d2f [t params]
      """Return Hessian of powerlaw"""
  (setv (, a b) (, (get params 0)
                   (get params 1)))
  (setv d2f (np.zeros (, (len params) (len params) (len t)))
        (get d2f (, 0 1)) (np.multiply (np.log t)
                                       (np.power t b))
        (get d2f (, 1 0)) (get d2f (, 0 1))
        (get d2f (, 1 0)) (np.multiply (* a (np.power (np.log t) 2))
                                       (np.power t b)))
  d2f)


;; string -> nil
(defn error [string &optional [stop True]]
  (sys.stderr.write (+ "- Error! " string))
  (if stop
      (sys.exit 2)))



;; string, string -> tuple of [np.array(1), np.array(2)]
(defn read-in-data [load-base-path]
  """Read in all M 2-column data files (trajectories) from path,
  return a list of N time points and an MxN matrix of corresponding
  trajectory values.

  (Code that is commented out allows for saving the data back as
  a single compressed binary file, useful if many load/read operations)
  """

  (setv trajectories (list)
        time (list))
  (unless load-base-path
      (error "Specify base path to folder with data to fit. Each file having two colums, <time, trajectory>.\n"))

  (setv file-names (glob.glob (+ load-base-path "/*")))

  (unless file-names
      (error "Specify <base path> to data files. No files found!\n"))

  (for [(, i file-name) (enumerate file-names)]
    (setv a-file (open file-name "r"))
    (.append trajectories (list))
    (for [line (.readlines a-file)]
      (setv li (.strip line))
      (when (not (.startswith li "#"))
        (setv tmp (list-comp
                    (float value)
                    (value (.split line))))
        (.append (get trajectories i) (get tmp 1))
        (if (= i 0)
            (.append time (get tmp 0)))))
    (.close a-file))
  (, (np.array time) (np.array trajectories)))

;; np.array(1), np.array(1) -> nil
(defn coff-det [y f]
  """Coefficient of determination, to check goodness-of-fit, from mean
  of data (MSD), and f(t, params).
  https://en.wikipedia.org/wiki/Coefficient_of_determination
  """
  (setv N (len y))

  ;; note, by y-mean here we mean the average over N, not M!
  (setv y-mean (/  (np.sum y) N)
        ;;SS-tot (np.sum (np.square (np.subtract y-mean y)))
        ;; using Hy-threadinig instead of above:
        SS-tot (-> (np.subtract y-mean y) np.square np.sum)
        SS-reg (-> (np.subtract y-mean f) np.square np.sum)
        SS-res (-> (np.subtract y      f) np.square np.sum))

  (print (+ "# Coefficient of determination:\t"
            (str (- 1 (/ SS-res SS-tot))))))


(defn main [args]
  "Perform a least squares including correlation in error on data in folder"

  (if (< (len args) 2)
      (do (print "usage: " (get args 0) " <data-path/>")
          (sys.exit 1)))

  ;; Read in data. Assumes path given is a folder where _all_ files
  ;; are 2 column files with x ("time") being column 1 and y
  ;; ("trajectory") being column 2..
  (setv (, time trajectories) (read-in-data (get args 1))
        (, M N) (np.shape trajectories))
  (setv min-method "nm")           ;; Use Nelder-Mead minimization method
  (setv guess (np.array '(.5 .5))) ;; Starting point in parameter space

  ;; OPTIONAL: Find index of first time to include in fitting procedure
  (setv starttime 200          ; first time point to include in fitting
        epsilon0 (- (get time 1) (get time 0)) ; time step size
        start-idx (int (/ starttime epsilon0)))

  ;;trajectories[:,start_idx:]
  (setv trajectories (get trajectories (, (slice None) (slice start-idx None))))
  (setv time (get time (, (slice start-idx None)))) ; does this work? should be [start-idx:]
  (setv (, M N) (np.shape trajectories))
  (assert (= N (len time)))
  (print "# Removed the first" start-idx "sampling times, new N =" N)
  ;; END

  ;; The analytical function to fit, its gradient, and hessian
  (.init wls-ice f df d2f)

  ;; Perform the actual fit
  (setv (, params sigma chi2-min) (wls-ice.fit time trajectories guess min-method))

  ;; RESULT:
  (print "# trajectories M=" M "\tsampling times N=" N "t_0=" (get time 0))
  (print "# Optimal param:\t" params)
  (print "# Sigma:\t" sigma)
  (print "# Chi-square value:\t" chi2-min)

  ;; Also get goodness-of-fitt parametes business
  (setv y-mean (wls-ice.compute-mean trajectories))
  (coff-det y-mean (f time params)))

(main sys.argv)
