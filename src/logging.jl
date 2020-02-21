using Logging
lg = ConsoleLogger(stderr,Logging.Debug)
global_logger(lg)


x = LandUse.runk(cpar = Dict(:S => 1.0, :L => 1.0, :kshare => [0.7,0.3], :K => 2, :θprop => [1.0,0.99]))
