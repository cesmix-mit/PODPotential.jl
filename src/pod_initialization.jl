using LinearAlgebra: eigen, Symmetric

# initialization only, can be done on CPU
function init2body(params)

    Φ = eigenvaluedecomposition!(2000,params)
    return Φ   
end

# initialization only, can be done on CPU
function eigenvaluedecomposition!(N,params)
    ns = params.counts.ns    
    rin = params.rin
    rcut = params.rcut

    # `snapshots` assumes S is a contiguous vector 
    S = Vector{Float64}(undef,N*ns)

    xij = [(rin + 1e-6) + (rcut - rin - 1e-6) * (i / (N - 1)) for i in 0:(N-1)]

    snapshots!(S,xij,N, params)
    
    # convert to a matrix for simplicity
    S = reshape(S,N,ns)

    # Compute S^TS (i.e. covariance matrix) 
    A = (1/N).*(S'*S)
    
    #Get eigenvectors
    _λ, Φ = eigen(Symmetric(A,:U), sortby=-)

    # project snapshots onto orthogonal basis
    Q = S*Φ
    
    # compute area to normalize the eignevectors
    # should be the same value at every entry?
    xij_diff  = [xij[i+1]-xij[i] for i in 1:N-1]
    for m in 1:ns
        area = 0.0
        for i in 1:N-1
            area += 0.5*xij_diff[i]*(Q[i,m]*Q[i,m] + Q[i+1,m]*Q[i+1,m])
        end
        Φ[:,m] ./= sqrt(area)
    end
   
    # enforce consistent signs for eigenvectors
    for m in 1:ns
        if Φ[m,m] < 0.0
            Φ[:,m] .*= -1
        end
    end

    return Φ
end

# initialization only, can be done on CPU
function snapshots!(S,xij,N, params)
    rcut = params.rcut 
    rin  = params.rin
    rmax = rcut - rin
    
    # βk values for k=1,2,3
    besselparams = [1e-3, 2.0, 4.0]

    Pα = params.besseldegree
    Pγ = params.inversedegree
    
    # iterate through the xij grid
    for n in eachindex(xij)
        dij = xij[n]

        r = dij - rin;

        #compute the normalized distance
        y = r/rmax;

        # cutoff function 
        y2 = y*y;
        y3 = 1.0 - y2*y;
        y4 = y3*y3 + 1e-6;
        y5 = sqrt(y4);
        y6 = exp(-1.0/y5);
        fcut = y6/exp(-1.0);
        
        # iterating over betas, fixed Pβ=3
        for k in 1:3 
            βk = besselparams[k]
            #this is useless?
            if(abs(βk) <= 1e-6)
                βk = 1e-3
            end

            x = (1.0 - exp(-βk*r/rmax))/(1.0-exp(-βk))

            for i in 1:Pα
                a = i*π
                b = (sqrt(2.0/(rmax))/i)

                nij = n + N*(i-1) + N*Pα*(k-1)

                S[nij] = b*fcut*sin(a*x)/r
            end
        end

        for i in 1:Pγ
            p = Pα*3 + (i-1)
            nij = n + N*p
            a = dij^i

            S[nij] = fcut/a
        end
    end
end    

function init3body(params::PODParams)
    nabf3 = params.counts.nabf3
    K3    = params.counts.K3
    P3    = params.P3

    pn3 = Vector{Int64}(undef,nabf3+1)
    pq3 = Vector{Int64}(undef,K3*2)
    pc3 = Vector{Int64}(undef,K3)

    init3bodyarray!(pn3, pq3, pc3, P3)
    return pn3, pq3, pc3
end

function init3bodyarray!(np, pq, pc, Pa)
    npa = (0,1,4,10,20,35,56,84,120,165,220,286,364,455)

    poly  = [
    (0, 0, 0, 0, 0, 1),
    (1, 0, 0, 1, 1, 1),
    (0, 1, 0, 1, 2, 1),
    (0, 0, 1, 1, 3, 1),
    (2, 0, 0, 2, 1, 1),
    (1, 1, 0, 2, 2, 2),
    (0, 2, 0, 3, 2, 1),
    (1, 0, 1, 2, 3, 2),
    (0, 1, 1, 3, 3, 2),
    (0, 0, 2, 4, 3, 1),
    (3, 0, 0, 5, 1, 1),
    (2, 1, 0, 5, 2, 3),
    (1, 2, 0, 6, 2, 3),
    (0, 3, 0, 7, 2, 1),
    (2, 0, 1, 5, 3, 3),
    (1, 1, 1, 6, 3, 6),
    (0, 2, 1, 7, 3, 3),
    (1, 0, 2, 8, 3, 3),
    (0, 1, 2, 9, 3, 3),
    (0, 0, 3, 10, 3, 1),
    (4, 0, 0, 11, 1, 1),
    (3, 1, 0, 11, 2, 4),
    (2, 2, 0, 12, 2, 6),
    (1, 3, 0, 13, 2, 4),
    (0, 4, 0, 14, 2, 1),
    (3, 0, 1, 11, 3, 4),
    (2, 1, 1, 12, 3, 12),
    (1, 2, 1, 13, 3, 12),
    (0, 3, 1, 14, 3, 4),
    (2, 0, 2, 15, 3, 6),
    (1, 1, 2, 16, 3, 12),
    (0, 2, 2, 17, 3, 6),
    (1, 0, 3, 18, 3, 4),
    (0, 1, 3, 19, 3, 4),
    (0, 0, 4, 20, 3, 1),
    (5, 0, 0, 21, 1, 1),
    (4, 1, 0, 21, 2, 5),
    (3, 2, 0, 22, 2, 10),
    (2, 3, 0, 23, 2, 10),
    (1, 4, 0, 24, 2, 5),
    (0, 5, 0, 25, 2, 1),
    (4, 0, 1, 21, 3, 5),
    (3, 1, 1, 22, 3, 20),
    (2, 2, 1, 23, 3, 30),
    (1, 3, 1, 24, 3, 20),
    (0, 4, 1, 25, 3, 5),
    (3, 0, 2, 26, 3, 10),
    (2, 1, 2, 27, 3, 30),
    (1, 2, 2, 28, 3, 30),
    (0, 3, 2, 29, 3, 10),
    (2, 0, 3, 30, 3, 10),
    (1, 1, 3, 31, 3, 20),
    (0, 2, 3, 32, 3, 10),
    (1, 0, 4, 33, 3, 5),
    (0, 1, 4, 34, 3, 5),
    (0, 0, 5, 35, 3, 1),
    (6, 0, 0, 36, 1, 1),
    (5, 1, 0, 36, 2, 6),
    (4, 2, 0, 37, 2, 15),
    (3, 3, 0, 38, 2, 20),
    (2, 4, 0, 39, 2, 15),
    (1, 5, 0, 40, 2, 6),
    (0, 6, 0, 41, 2, 1),
    (5, 0, 1, 36, 3, 6),
    (4, 1, 1, 37, 3, 30),
    (3, 2, 1, 38, 3, 60),
    (2, 3, 1, 39, 3, 60),
    (1, 4, 1, 40, 3, 30),
    (0, 5, 1, 41, 3, 6),
    (4, 0, 2, 42, 3, 15),
    (3, 1, 2, 43, 3, 60),
    (2, 2, 2, 44, 3, 90),
    (1, 3, 2, 45, 3, 60),
    (0, 4, 2, 46, 3, 15),
    (3, 0, 3, 47, 3, 20),
    (2, 1, 3, 48, 3, 60),
    (1, 2, 3, 49, 3, 60),
    (0, 3, 3, 50, 3, 20),
    (2, 0, 4, 51, 3, 15),
    (1, 1, 4, 52, 3, 30),
    (0, 2, 4, 53, 3, 15),
    (1, 0, 5, 54, 3, 6),
    (0, 1, 5, 55, 3, 6),
    (0, 0, 6, 56, 3, 1),
    (7, 0, 0, 57, 1, 1),
    (6, 1, 0, 57, 2, 7),
    (5, 2, 0, 58, 2, 21),
    (4, 3, 0, 59, 2, 35),
    (3, 4, 0, 60, 2, 35),
    (2, 5, 0, 61, 2, 21),
    (1, 6, 0, 62, 2, 7),
    (0, 7, 0, 63, 2, 1),
    (6, 0, 1, 57, 3, 7),
    (5, 1, 1, 58, 3, 42),
    (4, 2, 1, 59, 3, 105),
    (3, 3, 1, 60, 3, 140),
    (2, 4, 1, 61, 3, 105),
    (1, 5, 1, 62, 3, 42),
    (0, 6, 1, 63, 3, 7),
    (5, 0, 2, 64, 3, 21),
    (4, 1, 2, 65, 3, 105),
    (3, 2, 2, 66, 3, 210),
    (2, 3, 2, 67, 3, 210),
    (1, 4, 2, 68, 3, 105),
    (0, 5, 2, 69, 3, 21),
    (4, 0, 3, 70, 3, 35),
    (3, 1, 3, 71, 3, 140),
    (2, 2, 3, 72, 3, 210),
    (1, 3, 3, 73, 3, 140),
    (0, 4, 3, 74, 3, 35),
    (3, 0, 4, 75, 3, 35),
    (2, 1, 4, 76, 3, 105),
    (1, 2, 4, 77, 3, 105),
    (0, 3, 4, 78, 3, 35),
    (2, 0, 5, 79, 3, 21),
    (1, 1, 5, 80, 3, 42),
    (0, 2, 5, 81, 3, 21),
    (1, 0, 6, 82, 3, 7),
    (0, 1, 6, 83, 3, 7),
    (0, 0, 7, 84, 3, 1),
    (8, 0, 0, 85, 1, 1),
    (7, 1, 0, 85, 2, 8),
    (6, 2, 0, 86, 2, 28),
    (5, 3, 0, 87, 2, 56),
    (4, 4, 0, 88, 2, 70),
    (3, 5, 0, 89, 2, 56),
    (2, 6, 0, 90, 2, 28),
    (1, 7, 0, 91, 2, 8),
    (0, 8, 0, 92, 2, 1),
    (7, 0, 1, 85, 3, 8),
    (6, 1, 1, 86, 3, 56),
    (5, 2, 1, 87, 3, 168),
    (4, 3, 1, 88, 3, 280),
    (3, 4, 1, 89, 3, 280),
    (2, 5, 1, 90, 3, 168),
    (1, 6, 1, 91, 3, 56),
    (0, 7, 1, 92, 3, 8),
    (6, 0, 2, 93, 3, 28),
    (5, 1, 2, 94, 3, 168),
    (4, 2, 2, 95, 3, 420),
    (3, 3, 2, 96, 3, 560),
    (2, 4, 2, 97, 3, 420),
    (1, 5, 2, 98, 3, 168),
    (0, 6, 2, 99, 3, 28),
    (5, 0, 3, 100, 3, 56),
    (4, 1, 3, 101, 3, 280),
    (3, 2, 3, 102, 3, 560),
    (2, 3, 3, 103, 3, 560),
    (1, 4, 3, 104, 3, 280),
    (0, 5, 3, 105, 3, 56),
    (4, 0, 4, 106, 3, 70),
    (3, 1, 4, 107, 3, 280),
    (2, 2, 4, 108, 3, 420),
    (1, 3, 4, 109, 3, 280),
    (0, 4, 4, 110, 3, 70),
    (3, 0, 5, 111, 3, 56),
    (2, 1, 5, 112, 3, 168),
    (1, 2, 5, 113, 3, 168),
    (0, 3, 5, 114, 3, 56),
    (2, 0, 6, 115, 3, 28),
    (1, 1, 6, 116, 3, 56),
    (0, 2, 6, 117, 3, 28),
    (1, 0, 7, 118, 3, 8),
    (0, 1, 7, 119, 3, 8),
    (0, 0, 8, 120, 3, 1),
    (9, 0, 0, 121, 1, 1),
    (8, 1, 0, 121, 2, 9),
    (7, 2, 0, 122, 2, 36),
    (6, 3, 0, 123, 2, 84),
    (5, 4, 0, 124, 2, 126),
    (4, 5, 0, 125, 2, 126),
    (3, 6, 0, 126, 2, 84),
    (2, 7, 0, 127, 2, 36),
    (1, 8, 0, 128, 2, 9),
    (0, 9, 0, 129, 2, 1),
    (8, 0, 1, 121, 3, 9),
    (7, 1, 1, 122, 3, 72),
    (6, 2, 1, 123, 3, 252),
    (5, 3, 1, 124, 3, 504),
    (4, 4, 1, 125, 3, 630),
    (3, 5, 1, 126, 3, 504),
    (2, 6, 1, 127, 3, 252),
    (1, 7, 1, 128, 3, 72),
    (0, 8, 1, 129, 3, 9),
    (7, 0, 2, 130, 3, 36),
    (6, 1, 2, 131, 3, 252),
    (5, 2, 2, 132, 3, 756),
    (4, 3, 2, 133, 3, 1260),
    (3, 4, 2, 134, 3, 1260),
    (2, 5, 2, 135, 3, 756),
    (1, 6, 2, 136, 3, 252),
    (0, 7, 2, 137, 3, 36),
    (6, 0, 3, 138, 3, 84),
    (5, 1, 3, 139, 3, 504),
    (4, 2, 3, 140, 3, 1260),
    (3, 3, 3, 141, 3, 1680),
    (2, 4, 3, 142, 3, 1260),
    (1, 5, 3, 143, 3, 504),
    (0, 6, 3, 144, 3, 84),
    (5, 0, 4, 145, 3, 126),
    (4, 1, 4, 146, 3, 630),
    (3, 2, 4, 147, 3, 1260),
    (2, 3, 4, 148, 3, 1260),
    (1, 4, 4, 149, 3, 630),
    (0, 5, 4, 150, 3, 126),
    (4, 0, 5, 151, 3, 126),
    (3, 1, 5, 152, 3, 504),
    (2, 2, 5, 153, 3, 756),
    (1, 3, 5, 154, 3, 504),
    (0, 4, 5, 155, 3, 126),
    (3, 0, 6, 156, 3, 84),
    (2, 1, 6, 157, 3, 252),
    (1, 2, 6, 158, 3, 252),
    (0, 3, 6, 159, 3, 84),
    (2, 0, 7, 160, 3, 36),
    (1, 1, 7, 161, 3, 72),
    (0, 2, 7, 162, 3, 36),
    (1, 0, 8, 163, 3, 9),
    (0, 1, 8, 164, 3, 9),
    (0, 0, 9, 165, 3, 1),
    (10, 0, 0, 166, 1, 1),
    (9, 1, 0, 166, 2, 10),
    (8, 2, 0, 167, 2, 45),
    (7, 3, 0, 168, 2, 120),
    (6, 4, 0, 169, 2, 210),
    (5, 5, 0, 170, 2, 252),
    (4, 6, 0, 171, 2, 210),
    (3, 7, 0, 172, 2, 120),
    (2, 8, 0, 173, 2, 45),
    (1, 9, 0, 174, 2, 10),
    (0, 10, 0, 175, 2, 1),
    (9, 0, 1, 166, 3, 10),
    (8, 1, 1, 167, 3, 90),
    (7, 2, 1, 168, 3, 360),
    (6, 3, 1, 169, 3, 840),
    (5, 4, 1, 170, 3, 1260),
    (4, 5, 1, 171, 3, 1260),
    (3, 6, 1, 172, 3, 840),
    (2, 7, 1, 173, 3, 360),
    (1, 8, 1, 174, 3, 90),
    (0, 9, 1, 175, 3, 10),
    (8, 0, 2, 176, 3, 45),
    (7, 1, 2, 177, 3, 360),
    (6, 2, 2, 178, 3, 1260),
    (5, 3, 2, 179, 3, 2520),
    (4, 4, 2, 180, 3, 3150),
    (3, 5, 2, 181, 3, 2520),
    (2, 6, 2, 182, 3, 1260),
    (1, 7, 2, 183, 3, 360),
    (0, 8, 2, 184, 3, 45),
    (7, 0, 3, 185, 3, 120),
    (6, 1, 3, 186, 3, 840),
    (5, 2, 3, 187, 3, 2520),
    (4, 3, 3, 188, 3, 4200),
    (3, 4, 3, 189, 3, 4200),
    (2, 5, 3, 190, 3, 2520),
    (1, 6, 3, 191, 3, 840),
    (0, 7, 3, 192, 3, 120),
    (6, 0, 4, 193, 3, 210),
    (5, 1, 4, 194, 3, 1260),
    (4, 2, 4, 195, 3, 3150),
    (3, 3, 4, 196, 3, 4200),
    (2, 4, 4, 197, 3, 3150),
    (1, 5, 4, 198, 3, 1260),
    (0, 6, 4, 199, 3, 210),
    (5, 0, 5, 200, 3, 252),
    (4, 1, 5, 201, 3, 1260),
    (3, 2, 5, 202, 3, 2520),
    (2, 3, 5, 203, 3, 2520),
    (1, 4, 5, 204, 3, 1260),
    (0, 5, 5, 205, 3, 252),
    (4, 0, 6, 206, 3, 210),
    (3, 1, 6, 207, 3, 840),
    (2, 2, 6, 208, 3, 1260),
    (1, 3, 6, 209, 3, 840),
    (0, 4, 6, 210, 3, 210),
    (3, 0, 7, 211, 3, 120),
    (2, 1, 7, 212, 3, 360),
    (1, 2, 7, 213, 3, 360),
    (0, 3, 7, 214, 3, 120),
    (2, 0, 8, 215, 3, 45),
    (1, 1, 8, 216, 3, 90),
    (0, 2, 8, 217, 3, 45),
    (1, 0, 9, 218, 3, 10),
    (0, 1, 9, 219, 3, 10),
    (0, 0, 10, 220, 3, 1),
    (11, 0, 0, 221, 1, 1),
    (10, 1, 0, 221, 2, 11),
    (9, 2, 0, 222, 2, 55),
    (8, 3, 0, 223, 2, 165),
    (7, 4, 0, 224, 2, 330),
    (6, 5, 0, 225, 2, 462),
    (5, 6, 0, 226, 2, 462),
    (4, 7, 0, 227, 2, 330),
    (3, 8, 0, 228, 2, 165),
    (2, 9, 0, 229, 2, 55),
    (1, 10, 0, 230, 2, 11),
    (0, 11, 0, 231, 2, 1),
    (10, 0, 1, 221, 3, 11),
    (9, 1, 1, 222, 3, 110),
    (8, 2, 1, 223, 3, 495),
    (7, 3, 1, 224, 3, 1320),
    (6, 4, 1, 225, 3, 2310),
    (5, 5, 1, 226, 3, 2772),
    (4, 6, 1, 227, 3, 2310),
    (3, 7, 1, 228, 3, 1320),
    (2, 8, 1, 229, 3, 495),
    (1, 9, 1, 230, 3, 110),
    (0, 10, 1, 231, 3, 11),
    (9, 0, 2, 232, 3, 55),
    (8, 1, 2, 233, 3, 495),
    (7, 2, 2, 234, 3, 1980),
    (6, 3, 2, 235, 3, 4620),
    (5, 4, 2, 236, 3, 6930),
    (4, 5, 2, 237, 3, 6930),
    (3, 6, 2, 238, 3, 4620),
    (2, 7, 2, 239, 3, 1980),
    (1, 8, 2, 240, 3, 495),
    (0, 9, 2, 241, 3, 55),
    (8, 0, 3, 242, 3, 165),
    (7, 1, 3, 243, 3, 1320),
    (6, 2, 3, 244, 3, 4620),
    (5, 3, 3, 245, 3, 9240),
    (4, 4, 3, 246, 3, 11550),
    (3, 5, 3, 247, 3, 9240),
    (2, 6, 3, 248, 3, 4620),
    (1, 7, 3, 249, 3, 1320),
    (0, 8, 3, 250, 3, 165),
    (7, 0, 4, 251, 3, 330),
    (6, 1, 4, 252, 3, 2310),
    (5, 2, 4, 253, 3, 6930),
    (4, 3, 4, 254, 3, 11550),
    (3, 4, 4, 255, 3, 11550),
    (2, 5, 4, 256, 3, 6930),
    (1, 6, 4, 257, 3, 2310),
    (0, 7, 4, 258, 3, 330),
    (6, 0, 5, 259, 3, 462),
    (5, 1, 5, 260, 3, 2772),
    (4, 2, 5, 261, 3, 6930),
    (3, 3, 5, 262, 3, 9240),
    (2, 4, 5, 263, 3, 6930),
    (1, 5, 5, 264, 3, 2772),
    (0, 6, 5, 265, 3, 462),
    (5, 0, 6, 266, 3, 462),
    (4, 1, 6, 267, 3, 2310),
    (3, 2, 6, 268, 3, 4620),
    (2, 3, 6, 269, 3, 4620),
    (1, 4, 6, 270, 3, 2310),
    (0, 5, 6, 271, 3, 462),
    (4, 0, 7, 272, 3, 330),
    (3, 1, 7, 273, 3, 1320),
    (2, 2, 7, 274, 3, 1980),
    (1, 3, 7, 275, 3, 1320),
    (0, 4, 7, 276, 3, 330),
    (3, 0, 8, 277, 3, 165),
    (2, 1, 8, 278, 3, 495),
    (1, 2, 8, 279, 3, 495),
    (0, 3, 8, 280, 3, 165),
    (2, 0, 9, 281, 3, 55),
    (1, 1, 9, 282, 3, 110),
    (0, 2, 9, 283, 3, 55),
    (1, 0, 10, 284, 3, 11),
    (0, 1, 10, 285, 3, 11),
    (0, 0, 11, 286, 3, 1),
    (12, 0, 0, 287, 1, 1),
    (11, 1, 0, 287, 2, 12),
    (10, 2, 0, 288, 2, 66),
    (9, 3, 0, 289, 2, 220),
    (8, 4, 0, 290, 2, 495),
    (7, 5, 0, 291, 2, 792),
    (6, 6, 0, 292, 2, 924),
    (5, 7, 0, 293, 2, 792),
    (4, 8, 0, 294, 2, 495),
    (3, 9, 0, 295, 2, 220),
    (2, 10, 0, 296, 2, 66),
    (1, 11, 0, 297, 2, 12),
    (0, 12, 0, 298, 2, 1),
    (11, 0, 1, 287, 3, 12),
    (10, 1, 1, 288, 3, 132),
    (9, 2, 1, 289, 3, 660),
    (8, 3, 1, 290, 3, 1980),
    (7, 4, 1, 291, 3, 3960),
    (6, 5, 1, 292, 3, 5544),
    (5, 6, 1, 293, 3, 5544),
    (4, 7, 1, 294, 3, 3960),
    (3, 8, 1, 295, 3, 1980),
    (2, 9, 1, 296, 3, 660),
    (1, 10, 1, 297, 3, 132),
    (0, 11, 1, 298, 3, 12),
    (10, 0, 2, 299, 3, 66),
    (9, 1, 2, 300, 3, 660),
    (8, 2, 2, 301, 3, 2970),
    (7, 3, 2, 302, 3, 7920),
    (6, 4, 2, 303, 3, 13860),
    (5, 5, 2, 304, 3, 16632),
    (4, 6, 2, 305, 3, 13860),
    (3, 7, 2, 306, 3, 7920),
    (2, 8, 2, 307, 3, 2970),
    (1, 9, 2, 308, 3, 660),
    (0, 10, 2, 309, 3, 66),
    (9, 0, 3, 310, 3, 220),
    (8, 1, 3, 311, 3, 1980),
    (7, 2, 3, 312, 3, 7920),
    (6, 3, 3, 313, 3, 18480),
    (5, 4, 3, 314, 3, 27720),
    (4, 5, 3, 315, 3, 27720),
    (3, 6, 3, 316, 3, 18480),
    (2, 7, 3, 317, 3, 7920),
    (1, 8, 3, 318, 3, 1980),
    (0, 9, 3, 319, 3, 220),
    (8, 0, 4, 320, 3, 495),
    (7, 1, 4, 321, 3, 3960),
    (6, 2, 4, 322, 3, 13860),
    (5, 3, 4, 323, 3, 27720),
    (4, 4, 4, 324, 3, 34650),
    (3, 5, 4, 325, 3, 27720),
    (2, 6, 4, 326, 3, 13860),
    (1, 7, 4, 327, 3, 3960),
    (0, 8, 4, 328, 3, 495),
    (7, 0, 5, 329, 3, 792),
    (6, 1, 5, 330, 3, 5544),
    (5, 2, 5, 331, 3, 16632),
    (4, 3, 5, 332, 3, 27720),
    (3, 4, 5, 333, 3, 27720),
    (2, 5, 5, 334, 3, 16632),
    (1, 6, 5, 335, 3, 5544),
    (0, 7, 5, 336, 3, 792),
    (6, 0, 6, 337, 3, 924),
    (5, 1, 6, 338, 3, 5544),
    (4, 2, 6, 339, 3, 13860),
    (3, 3, 6, 340, 3, 18480),
    (2, 4, 6, 341, 3, 13860),
    (1, 5, 6, 342, 3, 5544),
    (0, 6, 6, 343, 3, 924),
    (5, 0, 7, 344, 3, 792),
    (4, 1, 7, 345, 3, 3960),
    (3, 2, 7, 346, 3, 7920),
    (2, 3, 7, 347, 3, 7920),
    (1, 4, 7, 348, 3, 3960),
    (0, 5, 7, 349, 3, 792),
    (4, 0, 8, 350, 3, 495),
    (3, 1, 8, 351, 3, 1980),
    (2, 2, 8, 352, 3, 2970),
    (1, 3, 8, 353, 3, 1980),
    (0, 4, 8, 354, 3, 495),
    (3, 0, 9, 355, 3, 220),
    (2, 1, 9, 356, 3, 660),
    (1, 2, 9, 357, 3, 660),
    (0, 3, 9, 358, 3, 220),
    (2, 0, 10, 359, 3, 66),
    (1, 1, 10, 360, 3, 132),
    (0, 2, 10, 361, 3, 66),
    (1, 0, 11, 362, 3, 12),
    (0, 1, 11, 363, 3, 12),
    (0, 0, 12, 364, 3, 1)]

    for i in 1:Pa+2
        np[i] = npa[i]
    end

    nmax = np[Pa+2]

    for i in 1:nmax
        pq[i]        = poly[i][4]
        pq[i+nmax]   = poly[i][5]
        pc[i]        = poly[i][6]
    end
end
