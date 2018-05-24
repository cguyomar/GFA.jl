function rc(seq::String)
    comp = Dict('A'=>'T', 'T'=>'A', 'G'=>'C', 'C'=>'G')
    return(reverse(map(x -> comp[x], seq)))
end
