for r0 = 0.1:0.1:1
    r = [r0:0.01:5];
   
    t = zeros(size(r,2),1);
    f = @(x)sqrt(x./((r0^2)*(x-r0)));
    for i = 1:1:size(r,2)
        t(i,:) =  integral(f,r0,r(i));
    end
    plot(r,t);
    hold on;
end