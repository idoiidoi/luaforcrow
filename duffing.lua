-- Parameters
local delta = 0.3
local gamma = 0.311
local omega = 1.1
local epsilon = 0.07 -- Coupling strength
local dt = 0.5 -- Increased time step size

-- State variables
local x1 = 0
local y1 = 0
local x2 = 0
local y2 = 0

function init()
    -- metro[1].time = dt
    metro[1].event = function(cycle)
        local t = cycle * dt
        local cos_omega_t = gamma * math.cos(omega * t)
        
        -- Calculate x1 and y1
        local k1_x1 = dt * y1
        local k1_y1 = dt * (x1 - x1^3 - delta*y1 + cos_omega_t - epsilon*(x1 - x2))
        local k2_x1 = dt * (y1 + k1_y1/2)
        local k2_y1 = dt * ((x1 + k1_x1/2) - (x1 + k1_x1/2)^3 - delta*(y1 + k1_y1/2) + cos_omega_t - epsilon*((x1 + k1_x1/2) - x2))
        local k3_x1 = dt * (y1 + k2_y1/2)
        local k3_y1 = dt * ((x1 + k2_x1/2) - (x1 + k2_x1/2)^3 - delta*(y1 + k2_y1/2) + cos_omega_t - epsilon*((x1 + k2_x1/2) - x2))
        local k4_x1 = dt * (y1 + k3_y1)
        local k4_y1 = dt * ((x1 + k3_x1) - (x1 + k3_x1)^3 - delta*(y1 + k3_y1) + cos_omega_t - epsilon*((x1 + k3_x1) - x2))
        
        x1 = x1 + 1/6 * (k1_x1 + 2*k2_x1 + 2*k3_x1 + k4_x1)
        y1 = y1 + 1/6 * (k1_y1 + 2*k2_y1 + 2*k3_y1 + k4_y1)
        
        -- Calculate x2 and y2
        local k1_x2 = dt * y2
        local k1_y2 = dt * (x2 - x2^3 - delta*y2 + cos_omega_t + epsilon*(x1 - x2))
        local k2_x2 = dt * (y2 + k1_y2/2)
        local k2_y2 = dt * ((x2 + k1_x2/2) - (x2 + k1_x2/2)^3 - delta*(y2 + k1_y2/2) + cos_omega_t + epsilon*(x1 - (x2 + k1_x2/2)))
        local k3_x2 = dt * (y2 + k2_y2/2)
        local k3_y2 = dt * ((x2 + k2_x2/2) - (x2 + k2_x2/2)^3 - delta*(y2 + k2_y2/2) + cos_omega_t + epsilon*(x1 - (x2 + k2_x2/2)))
        local k4_x2 = dt * (y2 + k3_y2)
        local k4_y2 = dt * ((x2 + k3_x2) - (x2 + k3_x2)^3 - delta*(y2 + k3_y2) + cos_omega_t + epsilon*(x1 - (x2 + k3_x2)))

        x2 = x2 + 1/6 * (k1_x2 + 2*k2_x2 + 2*k3_x2 + k4_x2)
        y2 = y2 + 1/6 * (k1_y2 + 2*k2_y2 + 2*k3_y2 + k4_y2)

        -- Output values
        output[1].slew = 0.1 -- Set the slew to 0 for no glide between voltages
        output[2].slew = 0.1
        metro[1].time = 0.01
        output[1].volts = 5*(x1+1)
        output[2].volts = 5*(y1+1)
        output[3].volts = 5*(x2+1)
        output[4].volts = 5*(y2+1)
    end

    metro[1]:start()
end
