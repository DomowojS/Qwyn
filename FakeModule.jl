module FakeModule
export FakeFunct

    function FakeFunct(P, x)
        P.x=x
        res = P.x[2]+P.y
        if res < 0
            return 10.0
        else
            return res
        end
    
    end

end