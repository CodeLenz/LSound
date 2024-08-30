
#
# Auxiliary functions defined in 
# https://julialinearalgebra.github.io/ArnoldiMethod.jl/stable/usage/02_spectral_transformations.html
#
struct ShiftAndInvert{TA,TB,TT}
    A_lu::TA
    B::TB
    temp::TT
end

function (M::ShiftAndInvert)(y,x)
    mul!(M.temp, M.B, x)
    ldiv!(y, M.A_lu, M.temp)
end

function construct_linear_map(A,B)
    a = ShiftAndInvert(lu(A),B,Vector{eltype(A)}(undef, size(A,1)))
    LinearMap{eltype(A)}(a, size(A,1), ismutating=true)
end

#
# Solve the generalized eigenvalue problem
#
"""
Evaluates the smaller nev eigenvalues for the generalized problem (A - Bλ)X = 0

    Solve_Eigen_(A::AbstractMatrix, B::AbstractMatrix, nev=4, positive=true, order=true, tol=1E-5, restarts=500)

A, B -> Abstract matrices
nev  -> Number of eigenvalues (2*nev are evaluated to account for possible negative values)   
positive -> only positive values are returned
order    -> eigenvalues are sorted in ascending order
tol and restarts are defined in ArnoldiMethod
tol_residue -> tolerance to check the norm of the residuer to EACH eingenpair
verbose -> if true, show some debug messages
ortho -> if true, check if eigenvectors are orthogonal

Returns
flag 
  1 if OK
 -1 if partialschur fails
 -2 if solution is not OK
 -3 if eigenvectors are not orthogonal (just if ortho==true)

eigenvalues (real vector)
eigenvectors (real matrix)

obs: If flag is -1 eigenvalues is zeros(1) and eigenvectors is zeros(1,1).

"""
function Solve_Eigen_(A::AbstractMatrix, B::AbstractMatrix, nev=4; positive=true, order=true, tol=1E-5, tol_residue=1E-4,
                      restarts=1000, verbose=false, ortho=false, cut_pos=0.0)

    # Return flag 
    #  1 if OK
    # -1 if partialschur fails
    # -2 if solution is not OK
    # -3 if eigenvectors are not orthogonal
  
    # Target the largest eigenvalues of the inverted problem
    decomp, history  = partialschur(construct_linear_map(A, B), nev=nev, tol=tol, restarts=restarts, which=:LM)

    # Test for convergency
    history.converged  || return -1, zeros(1), zeros(1,1) #throw("Solve_eigen:: $(string(history))")
     
    # Extract eigenvalues and eigenvectors
    λs_inv, X = partialeigen(decomp)

    # Eigenvalues have to be inverted to find the smallest eigenvalues of the non-inverted problem.
    λs = real.(1 ./ λs_inv)

    # There is no guarantee that eigenvalues are positive and are in crescent order. Thus, we 
    # can enforce this situation
    if positive
        
        # positive eigenvalues
        pp = λs.>cut_pos

        # Verify if we have at least one positive eigenvalues
        length(pp)>=1 || throw("Solve_Eigen_:: there are no positive eigenvalues - $(λs)")

        # Extract just the positive values
        λp = λs[pp] 

        # The same for the eigenvectors
        Xp = real.(X[:,pp])

    else
        λp = λs
        Xp = X
    end

    # Effective number of eigenvalues
    nλp = length(λp)

    if order
       p = sortperm(λp) 
    else
       p = collect(1:nλp)
    end
        
    # Alias to sorted solutions
    sorted_λp = @view λp[p]
    sorted_X  = real.(Xp[:,p])

    # Normalize sorted_X[:,i] - just to make it easier to verify the 
    # eigenvectors in the output of this function
    for i=1:nλp

        # Extreme values 
        mn = minimum(sorted_X[:,i])
        mx = maximum(sorted_X[:,i])

        # Take the absolute value of the most extreme
        normat = mx
        if abs(mn)>mx
            normat = mn
        end

        !isapprox(normat,0.0) || throw("Solve_Eigen_:: normalization factor is approx 0.0")
        
        sorted_X[:,i]  ./= normat
    end
    
    # It seems that everything is OK
    # Return both the eigenvalues and the eigenvectors
    return  1, sorted_λp, sorted_X

end
