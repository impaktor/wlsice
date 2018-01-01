(import sys scipy.optimize [numpy :as np])

(def f 0)
(def df 0)
(def d2f 0)

(defn error [string &optional [stop True]]
  "Auxiliary function for error printing"
  (sys.stderr.write (+ "-Error! " string "\n"))
  (if stop (sys.exit 2)))

;; function, function, function -> nil
(defn init [f_ df_ d2f_]
  "Set function to fit"
  (global f  df d2f)
  (def f f_)
  (def df df_)
  (def d2f d2f_))

;; np.array(2) -> np.array(1)
(defn compute-mean [trajectories]
  "Compute mean of trajectories"

  (setv (, M N) (np.shape trajectories)
        y-mean (np.zeros N))

  (for [i (range N)]
    (assoc y-mean i (/  (.sum np (get trajectories (, (slice None) i))) M)))
  y_mean)


;; np.array(2), np.array(1) -> np.array(2)
(defn make-Y [trajectories y-mean]
  """Take an M x N array of trajectories, and return an N x M matrix with:
  Y[n,m] = (y_n^(m) - \bar{y}_n)
  """
  (setv (, M N) (np.shape trajectories))

  ;; compute Y-matrix
  (setv Y (np.zeros (, N M)))
  (for [m (range M)]
    (setv (get Y (, (slice None) m))
           (np.subtract (get trajectories (, m (slice None))) y-mean)))
  Y)


;; np.array(2) -> np.array(2)
(defn make-covariance-from-Y [Y]
  "Take a matrix of dim N x M and compute the covariance matrix"
  (setv (, N M) (np.shape Y))
  (/ (np.dot Y (np.transpose Y))
     (- M 1.0)))

;; np.array(2) -> np.array(2), np.array(1)
(defn make-covariance-matrix [trajectories]
  "Make covariance matrix from an MxN array of trajectories"
  (setv (, M N) (np.shape trajectories)
        y-mean (compute-mean trajectories)
        C (make-covariance-from-Y (make-Y trajectories y-mean)))
  (, C y-mean))

;; np.array(2),  np.array(3), np.array(2), np.array(2), np.array(1) -> np.array(1)
(defn error-estimation [df d2f R C delta]
    """WLS-ICE error estimation, valid also for non-linear fitting. R could
    be R = inv(cov), or the diagonal of that, or some other symmetric
    matrix of our choosing. With N sampling points, and k parameters
    we have:

    df     is a  k x N      dimensional array
    d2f    is a  k x k x N  dimensional array
    delta  is a  N          dimensional array

    """
  (setv q (len df)
        N (len R)
        first-term (np.zeros (, q q)))

  (for [a (range q)
        b (range q)
        i (range N)
        j (range N)]
    (setv (get first-term (, a b))
          (+ (get first-term (, a b))
             (* 2
                (get d2f (, a b i))
                (get R (, i j))
                (get delta j)))))

  (setv second-term (* 2 (np.dot (np.dot df R)
                                 (np.transpose df))))
  (setv hessian (+ first-term
                   second-term))

  (setv H-inv (np.linalg.inv hessian)
        RCR (np.dot R (np.dot C R))
        error (np.zeros (, q q)))

  (for [a (range q)
        b (range q)
        c (range q)
        d (range q)]
    (setv df-RCR-df (np.dot (get df c)
                            (np.dot RCR
                                    (np.transpose (get df d)))))
    (setv (get error (, a b))
          (+ (get error (, a b))
             (* 4
                (get H-inv (, a c))
                df-RCR-df
                (get H-inv (, d b))))))
  (.diagonal error))


;; np.array(1), np.array(1), np.array(1), np.array(2) -> number
(defn chi2 [params t y R]
  """ Compute chi^2 function to be minimized, based on input
  (y-f)*R*(y-f)
  This gives us our fitted parameters in params.
  (If R=diag(inv(Cov)) we get least square fit)
  """
  (setv delta (np.subtract y (f t params)))
  (np.dot (np.dot delta R) delta))


;; np.array(1), np.array(1), np.array(1), np.array(2) -> np.array(1)
(defn chi2-jacobian [params t y R]
  """Jacobian of chi^2 function, i.e.
  -df/dparams[i]*[C^(-1)]*(y-f) + (y-f)*[C^(-1)]*(-df/dparams[i])
  return np.array of result for each i
  """
  (print "Params into jacobian:" params)
  (setv dY (df t params)
        Y (f t params)
        ans (np.zeros (len dY)))
  (for [i (range (len dY))]
    (setv (get ans i)
          (-  (* 2
                 (np.dot (np.dot (get dY i) R)
                         (np.subtract y Y))))))
  (print "Jacobian:\t%s\n" % ans)
  ans)

;; np.array(1), np.array(1), np.array(2), np.array(2), np,array(1), string,  ->
;;  tuple(np.array(1), np.array(1), number)
(defn minimize [t y R C guess-start min-method]
  """Minimize the chi-square function, to find optimal parameters and
    their esitmated error for function f"""
  (setv display False)

  (cond [(= min-method "bfgs")
         (setv result (scipy.optimize.minimize
                        chi2
                        guess-start
                        :args (, t y R)
                        :method "BFGS"
                        :jac chi2-jacobian
                        :options {"gtol" 1e-8 "disp" display}))]
        [(= min-method "nm")
         (setv result (scipy.optimize.minimize
                        chi2
                        guess-start
                        :args (, t y R)
                        :method "Nelder-Mead"
                        :options {"xtol" 1e-8 "disp" display}))]
        [True
         (error (+ "Unknown method: " min-method "\n"))])
  (setv params result.x
        chi-value (chi2 params t y R)
        sigma (np.sqrt (error-estimation (df t params) (d2f t params) R C
                                         (np.subtract (f t params) y))))
  (, params sigma chi-value))


;; np.array(1), np.array(2), string, np.array(1) -> tuple(np.array(1), np.array(1), number)
(defn fit [time trajectories guess &optional [min-method "nm"]]
  "Perform the correlated corrected least squares fit"

  (setv (, M N) (np.shape trajectories)
        (, C y-mean) (make-covariance-matrix trajectories)
        y-sigma (np.sqrt (np.diag C)))

  ;; LS-ICE
  (setv C (/ C M)
        R (np.linalg.inv (np.diag (np.diag C) 0)))

  ;; Do the parameter fitting, return: [opt_parameters, errors, chi2_min]
  (minimize time y-mean R C guess min-method))
