# an assitant function to find microstates such that i+j ≤ N-1
function Φ(N)
    X = Array{Int64,1}[]
    for i in 0:N-1
        for j in 0:N-1-i
            push!(X,[i,j])
        end
    end
    return X
end

# count the remaining number of bonds
function Nᵣ(macrostate)
    return ∑(macro_state2sizes(macrostate...))
end

function macro_state2sizes(i,j,k) # map the macrostates to its corresponding characterization of microstates
    nc = 0
    nl = 0
    nr = 0
    if i == 0
        # do nothing
    elseif i == 1
        nl += NL
    elseif i == 2
        nc += NL
    else
        println("i domain error")
    end
    if j == 0
        # do nothing
    elseif j == 1
        nr += NR
    elseif j == 2
        nc += NR
    else
        println("j domain error")
    end
    if k == 0
        # do nothing
    elseif k == 1
        nc += NM
    else
        println("k domain error")
    end
    list = [x for x in [nl,nc,nr] if x > 0]
    return list
end

function macro_state2exact_sizes(i,j,k) # map the macrostates to its corresponding characterization of microstates
    nc = 0
    nl = 0
    nr = 0
    if i == 0
        # do nothing
    elseif i == 1
        nl += NL
    elseif i == 2
        nl += NL
    else
        println("i domain error")
    end
    if j == 0
        # do nothing
    elseif j == 1
        nr += NR
    elseif j == 2
        nr += NR
    else
        println("j domain error")
    end
    if k == 0
        # do nothing
    elseif k == 1
        nc += NM
    else
        println("k domain error")
    end
    return [nl,nc,nr]
end

# map the unfacilitated state described by (X,Y,Z, microstates...)
# to the exact microstate of the form {0,1}^N
function state_mapping(state)
    X,Y,Z = state[1:3]
    sizes = macro_state2sizes(X,Y,Z)
    exact_sizes = macro_state2exact_sizes(X,Y,Z)
    if length(sizes) == 0
        exact_state = zeros(N)
    elseif length(sizes) == 1
        position_flag = [x > 0 for x in exact_sizes]
        undetermined_state = zeros(sizes[1])
        i,j = state[4]
        undetermined_state[i+1:sizes[1]-j] .= 1
        exact_state = Array{Int,1}()
        put_in_flag = true
        for i in 1:3
            flag = position_flag[i]
            if flag == 1 && put_in_flag
                exact_state = [exact_state;undetermined_state]
                put_in_flag = false
            end
            if flag == 0
                exact_state = [exact_state;zeros([NL,NM,NR][i])]
            end
        end
    elseif length(sizes) == 2
        undetermined_state_1 = zeros(sizes[1])
        i1,j1 = state[4]
        undetermined_state_1[i1+1:sizes[1]-j1] .= 1
        undetermined_state_2 = zeros(sizes[2])
        i2,j2 = state[5]
        undetermined_state_2[i2+1:sizes[2]-j2] .= 1
        undetermined_states = Array{Any,1}()
        push!(undetermined_states,undetermined_state_1)
        push!(undetermined_states,undetermined_state_2)
        exact_state = Array{Int,1}()
        cursor = 1
        position_flag = [x > 0 for x in exact_sizes]
        for i in 1:3
            flag = position_flag[i]
            if flag == 1
                if cursor < 3
                    exact_state = [exact_state; undetermined_states[cursor]]
                    cursor += 1
                end
            elseif flag == 0
                exact_state = [exact_state;zeros([NL,NM,NR][i])]
            end
        end
    elseif length(sizes) == 3
        L_state = zeros(NL)
        M_state = zeros(NM)
        R_state = zeros(NR)
        i1,j1 = state[4]
        i2,j2 = state[5]
        i3,j3 = state[6]
        L_state[i1+1:NL-j1] .= 1
        M_state[i2+1:NM-j2] .= 1
        R_state[i3+1:NR-j3] .= 1
        exact_state = [L_state; M_state; R_state]
    end
    return exact_state
end

############# Utility for visualizing the state space of facilitated case #############
# state_mapping_with_remodeller([1, 1, 0, [3, 0], [3, 0], [3, 0]]) = [
#    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
#    [3.0, 0.0]
#]
function state_mapping_with_remodeller(state)
    X,Y,Z = state[1:3]
    sizes = macro_state2sizes(X,Y,Z)
    exact_sizes = macro_state2exact_sizes(X,Y,Z)
    if length(sizes) == 0
        exact_state = zeros(N)
    elseif length(sizes) == 1
        position_flag = [x > 0 for x in exact_sizes]
        undetermined_state = zeros(sizes[1])
        i,j = state[4]
        undetermined_state[i+1:sizes[1]-j] .= 1
        exact_state = Array{Int,1}()
        put_in_flag = true
        for i in 1:3
            flag = position_flag[i]
            if flag == 1 && put_in_flag
                exact_state = [exact_state;undetermined_state]
                put_in_flag = false
            end
            if flag == 0
                exact_state = [exact_state;zeros([NL,NM,NR][i])]
            end
        end
    elseif length(sizes) == 2
        undetermined_state_1 = zeros(sizes[1])
        i1,j1 = state[4]
        undetermined_state_1[i1+1:sizes[1]-j1] .= 1
        undetermined_state_2 = zeros(sizes[2])
        i2,j2 = state[5]
        undetermined_state_2[i2+1:sizes[2]-j2] .= 1
        undetermined_states = Array{Any,1}()
        push!(undetermined_states,undetermined_state_1)
        push!(undetermined_states,undetermined_state_2)
        exact_state = Array{Int,1}()
        cursor = 1
        position_flag = [x > 0 for x in exact_sizes]
        for i in 1:3
            flag = position_flag[i]
            if flag == 1
                if cursor < 3
                    exact_state = [exact_state; undetermined_states[cursor]]
                    cursor += 1
                end
            elseif flag == 0
                exact_state = [exact_state;zeros([NL,NM,NR][i])]
            end
        end
    elseif length(sizes) == 3
        L_state = zeros(NL)
        M_state = zeros(NM)
        R_state = zeros(NR)
        i1,j1 = state[4]
        i2,j2 = state[5]
        i3,j3 = state[6]
        L_state[i1+1:NL-j1] .= 1
        M_state[i2+1:NM-j2] .= 1
        R_state[i3+1:NR-j3] .= 1
        exact_state = [L_state; M_state; R_state]
    end
    return [exact_state,state[end]]
end


function is_adjacent(x,y)
    X = S0[x]
    Y = S0[y]
    state_X = state_mapping(X)
    state_Y = state_mapping(Y)
    if ∑(abs.(state_X-state_Y)) > 1
        return false
    else
        return true
    end
end

# internal transition from y to x
function internal_on(x::Int,y::Int)
    X = S0[x]
    Y = S0[y]
    if length(X) == length(Y)
        Δ = X - Y
        if Δ[1:3] == [0,0,0]
            state_X = state_mapping(X)
            state_Y = state_mapping(Y)
            if ∑(abs.(state_X-state_Y)) == 1
                if ∑(state_X-state_Y) > 0
                    return 1
                else
                    return 0
                end
            else
                return 0
            end
        else
            return 0
        end
    else 
        return 0
    end
end

function internal_off(x::Int, y::Int)
    X = S0[x]
    Y = S0[y]
    if length(X) == length(Y)
        Δ = X - Y
        if Δ[1:3] == [0,0,0]
            state_X = state_mapping(X)
            state_Y = state_mapping(Y)
            if ∑(abs.(state_X-state_Y)) == 1
                if ∑(state_X-state_Y) > 0
                    return 0
                else
                    return 1
                end
            else
                return 0
            end
        else
            return 0
        end
    else 
        return 0
    end
end

# macrostates transition from y to x induced by disconnection of
# H2A-H2B dimer from (H3-H4)₂. It can split into one individual
# subunit and one two-component complex when splitted from a complete
# histone; or into two individual subunits.
function undocking(x::Int,y::Int,type::String)
    X = S0[x]
    Y = S0[y]
    state_X = state_mapping(X)
    state_Y = state_mapping(Y) # state of linkers between histone and DNA
    # for transitions induced by disconnection of H2A-H2B dimer from (H3-H4)₂, require state_X = state_Y.
    if state_X == state_Y # microstates have to be the same
        inds = [i for i in 1:3 if X[i] != Y[i]]
        ## individual subunits to undock and not leave
        if length(inds) == 1 # macrostates differ up to 1
            ind = inds[1]
            if X[ind] == 1 && Y[ind] == 2 && type == "simple"
                return 1
            elseif X[ind] == 0 && Y[ind] == 2 && type == "coupled"
                return 1
            else
                return 0
            end
        ## individual subunits to undock and leave
        elseif length(inds) == 2 && 3 ∈ inds
            ind = inds[1]
            if X[ind] == 1 && Y[ind] == 2 && X[3] == 0 && Y[3] == 1 && type == "coupled"# induced by disconnection of (H3-H4)₂ tetramer to H2A-H2B
                return 1
            else
                return 0
            end
        ## connected subunits (presumbly composed of two subunits) to undock and leave
        elseif length(inds) == 3 && sort(Y[1:3]) == [1,2,2] && sort(X[1:3]) == [0,0,1] && X[3] == 0 && type == "coupled"
            return 1
        else
            return 0
        end
    else
        return 0
    end
end


# rate of reverse reaction

function docking(x::Int,y::Int,type::String)
    X = S0[x]
    Y = S0[y]
    state_X = state_mapping(X)
    state_Y = state_mapping(Y) # state of linkers between histone and DNA
    # for transitions induced by connection of H2A-H2B dimer to (H3-H4)₂, require state_X = state_Y.
    if state_X == state_Y
        inds = [i for i in 1:3 if X[i] != Y[i]]
        if length(inds) == 1 # macrostates differ up to 1
            ind = inds[1]
            if X[ind] == 2 && Y[ind] == 1 && type == "simple" ## individual subunits to undock and not leave
                return 1
            elseif X[ind] == 2 && Y[ind] == 0  && type == "H2d"## individual subunits H2A-H2B to undock and leave 
                return 1
            else
                return 0
            end
        ## individual subunits to undock and leave
        elseif length(inds) == 2 && 3 ∈ inds
            ind = inds[1]
            if X[ind] == 2 && Y[ind] == 1 && X[3] == 1 && Y[3] == 0 && type == "H3_H4" # induced by disconnection of (H3-H4)₂ tetramer to H2A-H2B
                return 1 # may specify a different rate
            else
                return 0
            end
        ## connected subunits (presumbly composed of two subunits) to undock and leave
        elseif length(inds) == 3 && sort(X[1:3]) == [1,2,2] && sort(Y[1:3]) == [0,0,1] && Y[3] == 0 && type == "HEX"
            return 1
        else
            return 0
        end
    else
        return 0
    end
end

# if the connection between one subunits and other part of the histone is lost, and the connection between it and DNA has been lost previously, it may dissociate completely.
# this reaction is coupled to microstates change.
function dissociation(x::Int,y::Int)
    X = S0[x]
    Y = S0[y]
    state_X = state_mapping(X)
    state_Y = state_mapping(Y) # restriction on microstates is imposed by the fact that different macrostates have different restrictions on microstates.
    if ∑(abs.(state_X-state_Y)) == 1
        inds = [i for i in 1:3 if X[i] != Y[i]]
        if length(inds) == 1 && 3 ∉ inds # lost of one of (H2A-H2B)
            ind = inds[1]
            if X[ind] == 0 && Y[ind] == 1
                return 1
            else
                return 0
            end
        elseif length(inds) == 1 && 3 ∈ inds # lost of (H3-H4)₂, note that it CAN occur when X = [0,0,0], Y = [0,1,0] and CANNOT occur when either one of H2A-H2B attaches to (H3-H4)₂.
            ind = inds[1]
            if maximum(X[1:3]) < 2 && maximum(Y[1:3]) < 2 && X[3] == 0 && Y[3] == 1 # no connection is established between (H3-H4)₂ and (H2A-H2B)
                return 1
            else
                return 0
            end
        elseif length(inds) == 2 && sort(Y[inds]) == [1,2] && sort(X[inds]) == [0,0] && 3 ∈ inds 
            return 1
        elseif length(inds) == 3 && Y[inds] == [2,2,1] && X[inds] == [0,0,0]
            return 1
        else
            return 0
        end
    else
        return 0
    end
end

# reverse reaction of dissociation
function association(x::Int,y::Int,type::String)
    X = S0[x]
    Y = S0[y]
    state_X = state_mapping(X)
    state_Y = state_mapping(Y) # restriction on microstates is imposed by the fact that different macrostates have different restrictions on microstates.
    if ∑(abs.(state_X-state_Y)) == 1
        inds = [i for i in 1:3 if X[i] != Y[i]]
        if length(inds) == 1 && 3 ∉ inds
            ind = inds[1]
            if X[ind] == 1 && Y[ind] == 0 && type == "H2d"
                return 1
            else
                return 0
            end
        elseif length(inds) == 1 && 3 ∈ inds # lost of (H3-H4)₂
            ind = inds[1]
            if maximum(X[1:3]) < 2 && maximum(Y[1:3]) < 2 && X[3] == 1 && Y[3] == 0 && type == "H3_H4" # no connection is established between (H3-H4)₂ and (H2A-H2B)
                return 1
            else
                return 0
            end
        elseif length(inds) == 2 && sort(X[inds]) == [1,2] && sort(Y[inds]) == [0,0] && 3 ∈ inds && type == "HEX" # (H3-H4)₂-H2A-H2B hexamer dissociate/associate
            return 1
        elseif length(inds) == 3 && X[1:3]== [2,2,1] && Y[1:3] == [0,0,0] && type == "OCT"
            return 1
        else
            return 0
        end
    else
        return 0
    end
end

function is_adjacent_with_remodeller(x,y)
    X = S0[x]
    Y = S0[y]
    state_X, helicase_X = state_mapping_with_remodeller(X)
    state_Y, helicase_Y  = state_mapping_with_remodeller(Y)
    if ∑(abs.(state_X-state_Y))+∑(abs.(helicase_X-helicase_Y)) > 1
        return false
    else
        return true
    end
end

# internal transition from y to x
function internal_on_with_remodeller(x::Int,y::Int)
    X = S0[x]
    Y = S0[y]
    if length(X) == length(Y)
        Δ = X - Y
        if Δ[1:3] == [0,0,0]
        state_X, helicase_X = state_mapping_with_remodeller(X)
        state_Y, helicase_Y  = state_mapping_with_remodeller(Y)
            if ∑(abs.(state_X-state_Y))+∑(abs.(helicase_X-helicase_Y)) == 1
                if ∑(state_X-state_Y) > 0
                    return 1
                else
                    return 0
                end
            else
                return 0
            end
        else
            return 0
        end
    else 
        return 0
    end
end

function internal_off_with_remodeller(x::Int, y::Int)
    X = S0[x]
    Y = S0[y]
    if length(X) == length(Y)
        Δ = X - Y
        if Δ[1:3] == [0,0,0]
            state_X, helicase_X = state_mapping_with_remodeller(X)
            state_Y, helicase_Y  = state_mapping_with_remodeller(Y)
            if ∑(abs.(state_X-state_Y))+∑(abs.(helicase_X-helicase_Y)) == 1
                if ∑(state_X-state_Y) < 0
                    return 1
                else
                    return 0
                end
            else
                return 0
            end
        else
            return 0
        end
    else 
        return 0
    end
end

function helicase_on(x::Int,y::Int)
    X = S0[x]
    Y = S0[y]
    if length(X) == length(Y)
        Δ = X - Y
        if Δ[1:3] == [0,0,0]
        state_X, helicase_X = state_mapping_with_remodeller(X)
        state_Y, helicase_Y  = state_mapping_with_remodeller(Y)
            if ∑(abs.(state_X-state_Y))+∑(abs.(helicase_X-helicase_Y)) == 1
                if ∑(helicase_X-helicase_Y) > 0
                    return 1
                else
                    return 0
                end
            else
                return 0
            end
        else
            return 0
        end
    else 
        return 0
    end
end

function helicase_off(x::Int,y::Int)
    X = S0[x]
    Y = S0[y]
    if length(X) == length(Y)
        Δ = X - Y
        if Δ[1:3] == [0,0,0]
        state_X, helicase_X = state_mapping_with_remodeller(X)
        state_Y, helicase_Y  = state_mapping_with_remodeller(Y)
            if ∑(abs.(state_X-state_Y))+∑(abs.(helicase_X-helicase_Y)) == 1
                if ∑(helicase_X-helicase_Y) < 0
                    return 1
                else
                    return 0
                end
            else
                return 0
            end
        else
            return 0
        end
    else 
        return 0
    end
end

# macrostates transition from y to x induced by disconnection of
# H2A-H2B dimer from (H3-H4)₂. It can split into one individual
# subunit and one two-component complex when splitted from a complete
# histone; or into two individual subunits.
function undocking_with_remodeller(x::Int,y::Int,type::String)
    X = S0[x]
    Y = S0[y]
    state_X, helicase_X = state_mapping_with_remodeller(X)
    state_Y, helicase_Y  = state_mapping_with_remodeller(Y) # state of linkers between histone and DNA
    # for transitions induced by disconnection of H2A-H2B dimer from (H3-H4)₂, require state_X = state_Y.
    if state_X == state_Y && helicase_X == helicase_Y # microstates have to be the same
        inds = [i for i in 1:3 if X[i] != Y[i]]
        ## individual subunits to undock and not leave
        if length(inds) == 1 # macrostates differ up to 1
            ind = inds[1]
            if X[ind] == 1 && Y[ind] == 2 && type == "simple"
                return 1
            elseif X[ind] == 0 && Y[ind] == 2 && type == "coupled"
                return 1
            else
                return 0
            end
        ## individual subunits to undock and leave
        elseif length(inds) == 2 && 3 ∈ inds
            ind = inds[1]
            if X[ind] == 1 && Y[ind] == 2 && X[3] == 0 && Y[3] == 1 && type == "coupled"# induced by disconnection of (H3-H4)₂ tetramer to H2A-H2B
                return 1
            else
                return 0
            end
        ## connected subunits (presumbly composed of two subunits) to undock and leave
        elseif length(inds) == 3 && sort(Y[1:3]) == [1,2,2] && sort(X[1:3]) == [0,0,1] && X[3] == 0 && type == "coupled"
            return 1
        else
            return 0
        end
    else
        return 0
    end
end


# rate of reverse reaction

function docking_with_remodeller(x::Int,y::Int,type::String)
    X = S0[x]
    Y = S0[y]
    state_X, helicase_X = state_mapping_with_remodeller(X)
    state_Y, helicase_Y  = state_mapping_with_remodeller(Y) # state of linkers between histone and DNA
    # for transitions induced by connection of H2A-H2B dimer to (H3-H4)₂, require state_X = state_Y.
    if state_X == state_Y  && helicase_X == helicase_Y 
        inds = [i for i in 1:3 if X[i] != Y[i]]
        if length(inds) == 1 # macrostates differ up to 1
            ind = inds[1]
            if X[ind] == 2 && Y[ind] == 1 && type == "simple" ## individual subunits to undock and not leave
                return 1
            elseif X[ind] == 2 && Y[ind] == 0  && type == "H2d"## individual subunits H2A-H2B to undock and leave 
                return 1
            else
                return 0
            end
        ## individual subunits to undock and leave
        elseif length(inds) == 2 && 3 ∈ inds
            ind = inds[1]
            if X[ind] == 2 && Y[ind] == 1 && X[3] == 1 && Y[3] == 0 && type == "H3_H4" # induced by disconnection of (H3-H4)₂ tetramer to H2A-H2B
                return 1 # may specify a different rate
            else
                return 0
            end
        ## connected subunits (presumbly composed of two subunits) to undock and leave
        elseif length(inds) == 3 && sort(X[1:3]) == [1,2,2] && sort(Y[1:3]) == [0,0,1] && Y[3] == 0 && type == "HEX"
            return 1
        else
            return 0
        end
    else
        return 0
    end
end

# if the connection between one subunits and other part of the histone is lost, and the connection between it and DNA has been lost previously, it may dissociate completely.
# this reaction is coupled to microstates change.
function dissociation_with_remodeller_competing(x::Int, y::Int)
    # used only for evaluating the competing pathways
    # designed to retain history of histone state upon full detachment
    X = S0[x]
    Y = S0[y]
    state_X, helicase_X = state_mapping_with_remodeller(X)
    state_Y, helicase_Y  = state_mapping_with_remodeller(Y) # restriction on microstates is imposed by the fact that different macrostates have different restrictions on microstates.
    # special treatment for the full detachment of histone event
    # this retains the history of the histone state before detaching
    if ∑(abs.(state_X)) == 0 && ∑(abs.(state_Y)) == 1
        if X[1:3] == Y[1:3]
            return 1
        else
            return 0
        end
    end
    if ∑(abs.(state_X-state_Y)) == 1  && helicase_X == helicase_Y 
        inds = [i for i in 1:3 if X[i] != Y[i]]
        if length(inds) == 1 && 3 ∉ inds # lost of one of (H2A-H2B)
            ind = inds[1]
            if X[ind] == 0 && Y[ind] == 1
                return 1
            else
                return 0
            end
        elseif length(inds) == 1 && 3 ∈ inds # lost of (H3-H4)₂, note that it CAN occur when X = [0,0,0], Y = [0,1,0] and CANNOT occur when either one of H2A-H2B attaches to (H3-H4)₂.
            ind = inds[1]
            if maximum(X[1:3]) < 2 && maximum(Y[1:3]) < 2 && X[3] == 0 && Y[3] == 1 # no connection is established between (H3-H4)₂ and (H2A-H2B)
                return 1
            else
                return 0
            end
        elseif length(inds) == 2 && sort(Y[inds]) == [1,2] && sort(X[inds]) == [0,0] && 3 ∈ inds 
            return 1
        elseif length(inds) == 3 && Y[inds] == [2,2,1] && X[inds] == [0,0,0]
            return 1
        else
            return 0
        end
    else
        return 0
    end
end


function dissociation_with_remodeller(x::Int,y::Int)
    X = S0[x]
    Y = S0[y]
    state_X, helicase_X = state_mapping_with_remodeller(X)
    state_Y, helicase_Y  = state_mapping_with_remodeller(Y) # restriction on microstates is imposed by the fact that different macrostates have different restrictions on microstates.
    if ∑(abs.(state_X-state_Y)) == 1  && helicase_X == helicase_Y 
        inds = [i for i in 1:3 if X[i] != Y[i]]
        if length(inds) == 1 && 3 ∉ inds # lost of one of (H2A-H2B)
            ind = inds[1]
            if X[ind] == 0 && Y[ind] == 1
                return 1
            else
                return 0
            end
        elseif length(inds) == 1 && 3 ∈ inds # lost of (H3-H4)₂, note that it CAN occur when X = [0,0,0], Y = [0,1,0] and CANNOT occur when either one of H2A-H2B attaches to (H3-H4)₂.
            ind = inds[1]
            if maximum(X[1:3]) < 2 && maximum(Y[1:3]) < 2 && X[3] == 0 && Y[3] == 1 # no connection is established between (H3-H4)₂ and (H2A-H2B)
                return 1
            else
                return 0
            end
        elseif length(inds) == 2 && sort(Y[inds]) == [1,2] && sort(X[inds]) == [0,0] && 3 ∈ inds 
            return 1
        elseif length(inds) == 3 && Y[inds] == [2,2,1] && X[inds] == [0,0,0]
            return 1
        else
            return 0
        end
    else
        return 0
    end
end

# reverse reaction of dissociation_with_remodeller
function association_with_remodeller(x::Int,y::Int,type::String)
    X = S0[x] ## X is the state of the histone *after* the reaction
    Y = S0[y] ## Y is the state of the histone *before* the reaction
    state_X, helicase_X = state_mapping_with_remodeller(X)
    state_Y, helicase_Y  = state_mapping_with_remodeller(Y) # restriction on microstates is imposed by the fact that different macrostates have different restrictions on microstates.
    if ∑(abs.(state_X-state_Y)) == 1  && helicase_X == helicase_Y 
        inds = [i for i in 1:3 if X[i] != Y[i]]
        if length(inds) == 1 && 3 ∉ inds
            ind = inds[1]
            if X[ind] == 1 && Y[ind] == 0 && type == "H2d"
                return 1
            else
                return 0
            end
        elseif length(inds) == 1 && 3 ∈ inds # lost of (H3-H4)₂
            ind = inds[1]
            if maximum(X[1:3]) < 2 && maximum(Y[1:3]) < 2 && X[3] == 1 && Y[3] == 0 && type == "H3_H4" # no connection is established between (H3-H4)₂ and (H2A-H2B)
                return 1
            else
                return 0
            end
        elseif length(inds) == 2 && sort(X[inds]) == [1,2] && sort(Y[inds]) == [0,0] && 3 ∈ inds && type == "HEX" # (H3-H4)₂-H2A-H2B hexamer dissociate/associate
            return 1
        elseif length(inds) == 3 && X[1:3]== [2,2,1] && Y[1:3] == [0,0,0] && type == "OCT"
            return 1
        else
            return 0
        end
    else
        return 0
    end
end


# another utility function for the construction of the state space
# in this model, we assume that the remodeler can only attack from ends toward the middle
# so we just need only 1 pair of indices to specify the position of the remodeler 
# the behavior of the remodeler is best described by *helicase*.
function get_plausible_helicase_postion(microstate)
    list = Array{Int,1}[]
    i=j= 0
    while i+1 < N && microstate[i+1] == 0
        i += 1
    end
    # restrict to the case where at least one site is occupied
    while j+1 < N && microstate[N-j] == 0
        j += 1
    end
    for a in 0:1:i, b in 0:1:j
        if a + b < N # only meant to deal with some exceptions, like microstate == 𝟎
            push!(list,[a,b])
        end
    end
    return list
end
