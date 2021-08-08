#
# Types for orbital trajectories
#

const Trajectory =  DifferentialEquations.SciMLBase.ODESolution{T,N,uType,uType2,DType,tType,rateType,P,A,IType,DE} where {
    T, N, uType <: Vector{<:States.AbstractState}, uType2,
    DType, tType, rateType, P, A, IType, DE 
}