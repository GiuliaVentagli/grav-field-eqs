using Symbolics, LinearAlgebra

#Covariant partial derivative of covariant metric
function derMet(Met::AbstractMatrix{T}) where T
    Dg = Matrix{Matrix}(undef,4,4) 
    for i = 1:n
        for j = 1:n
            Inset = Matrix{Num}(undef,1,4) 
            for k = 1:n
                Inset[k] = transpose(Symbolics.gradient(Met[i,j],[t,r,θ,ϕ]))[k]
            end
            Dg[i,j] = Inset
        end
    end
    return Dg
end

#Christoffel symbol udd
function christoffel(Met::AbstractMatrix{T}, Metuu::AbstractMatrix{T}) where T
    chris = Matrix{Matrix}(undef,4,4)
    for i = 1:n #u
        for k = 1:n #d
            Inset = Matrix{Num}(undef,1,4)
            for l = 1:n #d
                temp = 0
                for j = 1:n
                    temp += 1/2*Metuu[i,j] * (derMet(Met)[j,k][l]+derMet(Met)[j,l][k]-
                        derMet(Met)[k,l][j])
                end
                Inset[l] = temp
            end
            chris[i,k] = Inset
        end
    end
    return chris
end

#covariant partial derivative of christoffel symbol
function derChris(Met::AbstractMatrix{T}) where T
    DΓ = Matrix{Matrix}(undef,4,4) 
    for i = 1:n
        for j = 1:n
            Inset = Matrix{Num}(undef,4,4)
            for l = 1:n
                for k = 1:n
                    Inset[l,k] = transpose(
                            Symbolics.gradient(Met[i,j][l],[t,r,θ,ϕ]))[k]
                end
            end
            DΓ[i,j] = Inset
        end
    end
    return DΓ
end

#riemann tensor uddd
function riemann(chris::AbstractMatrix{Matrix})
    riem = Matrix{Matrix}(undef,4,4)
    for i = 1:n #u
        for j = 1:n #d
            Inset = Matrix{Num}(undef,4,4)
            for k = 1:n #d
                for l = 1:n #d
                    temp = derChris(chris)[i,j][l,k] -
                                derChris(chris)[i,j][k,l]
                    for s = 1:n
                        temp += chris[s,j][l] * chris[i,k][s] -
                                chris[s,j][k] * chris[i,l][s]
                    end
                    Inset[k,l] = temp
                end
            end
            riem[i,j] = Inset
        end
    end
    return riem
end

#Ricci tensor dd
function ricciTens(riem::AbstractMatrix{Matrix})
    ricciT = Matrix{Num}(undef,4,4)
    for j = 1:n #d
        for l = 1:n #d
            temp = 0
            for i = 1:n
                temp += riem[i,j][i,l]
            end
            ricciT[j,l] = temp
        end
    end
    return ricciT
end

#Ricci scalar
function ricciScal(Metuu::AbstractMatrix{T},ricciT::AbstractMatrix{T}) where T
    ricci = 0
    for i = 1:n
        for j = 1:n
            ricci += Metuu[i,j] * ricciT[i,j]
        end
    end
    return simplify(ricci; expand=true)
end

#Gauss-Bonnet invariant
function gaussBonnet(Met::AbstractMatrix{T},Metuu::AbstractMatrix{T},riem::AbstractMatrix{Matrix},
        ricciT::AbstractMatrix{T},scalar::T) where T

    gaussBonnet = scalar^2

    for i = 1:n
        for j = 1:n
            for k = 1:n
                for l = 1:n
                    gaussBonnet += -4 * Metuu[i,j] * Metuu[k,l] * ricciT[i,k] * ricciT[j,l] 
                end
            end
        end
    end

    for m = 1:n
        for j = 1:n
            for k = 1:n
                for l = 1:n
                    for i = 1:n
                        for x = 1:n
                            for c =1:n
                                for d = 1:n
                                    gaussBonnet += Met[m,i] * Metuu[j,x] * Metuu[k,c] * 
                                                Metuu[l,d] * riem[i,j][k,l] * riem[m,x][c,d] 
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return gaussBonnet
end

#stress-energy tensor T^{ab}
function stressEnergy(Metuu::AbstractMatrix{T}, 
        vel::AbstractMatrix{T},ϵ,p) where T
    stressEn = Matrix{Num}(undef,4,4)
    for i = 1:n
        for j = 1:n
            stressEn[i,j] = (ϵ + p) * vel[i] * vel[j] + p * Metuu[i,j]
        end
    end
    return stressEn
end

#stress-energy tensor T_{ab}
function stressEnergydd(Met::AbstractMatrix{T}, Tuu::AbstractMatrix{T}) where T
    stressEndd = Matrix{Num}(undef,4,4)
    for i = 1:n
        for j = 1:n
            temp = 0
            for k = 1:n
                for l = 1:n
                    temp += Met[i,k] * Met[j,l] * Tuu[k,l]
                end
            end
            stressEndd[i,j] = temp
        end
    end
    return stressEndd
end

#Conservation of stress-energy tensor
function DerStressEn(chris::AbstractMatrix{Matrix}, 
        Tuu::AbstractMatrix{T}) where T
    DStressEn = Matrix{Num}(undef,1,4)
    for j = 1:n
        temp = 0
        for i = 1:n 
            temp += transpose(Symbolics.gradient(Tuu[i,j],[t,r,θ,ϕ]))[i]
            for k = 1:n
                temp += chris[i,k][i] * Tuu[k,j] + chris[j,i][k] * Tuu[i,k]
            end
        end
        DStressEn[j] = temp
    end
    return DStressEn
end

#D_a@D_b@φ where D is the covariant derivative and φ is a scalar field
function derDerScal(chris::AbstractMatrix{Matrix}, scal::T) where T
    DDscal = Matrix{Num}(undef,4,4)
    for i = 1:n
        for j = 1:n
            DDscal[i,j] = 
                    Symbolics.jacobian(Symbolics.jacobian([scal],[t r θ ϕ]),[t r θ ϕ])[i,j]
            for k = 1:n
                DDscal[i,j] += -chris[k,i][j] * (Symbolics.jacobian([scal],[t r θ ϕ])[k])
            end
        end
    end
    return DDscal
end

#g^{ab}*D_a@D_b@φ where D is the covariant derivative and φ is a scalar field
function boxScal(Metuu::AbstractMatrix{T}, chris::AbstractMatrix{Matrix},
        scal::T) where T
    DDscal = 0
    for i = 1:n
        for j = 1:n
            DDscal += Metuu[i,j] * 
                        Symbolics.jacobian(Symbolics.jacobian([scal],[t r θ ϕ]),[t r θ ϕ])[i,j]
            for k = 1:n
                DDscal += - Metuu[i,j] * 
                            chris[k,i][j] * (Symbolics.jacobian([scal],[t r θ ϕ])[k])
            end
        end
    end
    return simplify(DDscal)
end;

#D_a@φ*D_b@φ where D is the covariant derivative and φ is a scalar field
function dScaldScal(scal::T) where T
    DscalDscal = Matrix{Num}(undef,4,4)
    for i = 1:n
        for j = 1:n
            DscalDscal[i,j] = Symbolics.jacobian([scal],[t r θ ϕ])[i] * 
                                Symbolics.jacobian([scal],[t r θ ϕ])[j]
        end
    end 
    return DscalDscal
end

#g^{ab}*D_a@φ*D_b@φ where D is the covariant derivative and φ is a scalar field
function dScaldScalContr(Metuu::AbstractMatrix{T}, scal::T) where T
    DscalDscal = 0
    for i = 1:n
        for j = 1:n
            DscalDscal += Metuu[i,j] * Symbolics.jacobian([scal],[t r θ ϕ])[i] * 
                            Symbolics.jacobian([scal],[t r θ ϕ])[j]
        end
    end    
    return DscalDscal
end;