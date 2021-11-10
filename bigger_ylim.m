function bigger_ylim 

    SF = 0.1; 
    ylimits = [get(gca, 'ylim')]; 
    yrange = ylimits(2) - ylimits(1); 
    new_ylim = [ ylimits(1) - SF*yrange, ylimits(2) + SF*yrange ]; 
    set(gca, 'ylim', new_ylim); 

end 