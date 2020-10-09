using NLsolve
using Makie
AbstractPlotting.__init__()
n = 12
function on_curve(n,d,f)
    x = zeros(n)
    y = zeros(n)
    y[1] = f(x[1])
    for i = 2:n
        function func!(func,unk)
            xi,yi = unk
            dy = yi - y[i-1]
            dx = xi - x[i-1]
            func[1] = dy^2 + dx^2 - d
            func[2] = yi-f(xi)
        end
        soli = nlsolve(func!,[x[i-1]+d,f(x[i-1]+d)])
        x[i],y[i] = soli.zero
    end
    x,y
end
x,y = on_curve(12,0.1,sin)
plot(x,y)
