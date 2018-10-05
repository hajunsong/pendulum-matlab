function define_Y_vector

    global body num_body
    global Y
    
    if num_body == 1
        Y(1,1) = body.qi;
        Y(2,1) = body.dqi;
    else
        for i = 1 : num_body
            Y(i,1) = body(i).qi;
        end
        for i = 1 : num_body
            Y(i+num_body,1) = body(i).dqi;
        end
    end

end