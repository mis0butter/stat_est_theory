function bigger_xlim 

    SF = 0.1; 
    xlimits = [get(gca, 'xlim')]; 
    xrange = xlimits(2) - xlimits(1); 
    new_xlim = [ xlimits(1) - SF*xrange, xlimits(2) + SF*xrange ]; 
    set(gca, 'xlim', new_xlim); 

end 