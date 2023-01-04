function f = timeObjFunc(x, h)
% TIMEOBJFUNC objective function
f=0;
for i=1:length(x)-1  
    f = f + 1/(sqrt(x(i)) + sqrt(x(i+1)));
end
f = 2*h*f;