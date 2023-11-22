using FFTW,Statistics

"""
    noise.power_law_noise(length::Integer, NoiseType)

Generate a series of normalized power-law noise  

```math
    s(f) = h * f^(NoiseType) 
```
NoiseType = -1 (flicker) by default

"""
power_law_noise(len::Integer) = power_law_noise(len, -1)
power_law_noise(len::Integer, NoiseType) = power_law_noise(len, NoiseType)

function power_law_noise(len::Integer, NoiseType)
    white_noise = randn(len)
    white_noise_fft = FFTW.rfft(white_noise)
    frequencies = FFTW.rfftfreq(len)
    # Adjust the amplitude of each frequency component according to the power law
    # Avoid division by zero for the first element
    power_law_spectrum = white_noise_fft ./ (frequencies .+ 1e-10).^ (NoiseType / 2)
    # Inverse FFT to get noise in time domain
    power_law_noise = real(FFTW.irfft(power_law_spectrum, len))
    # Center and normalize the noise
    power_law_noise = (power_law_noise .- mean(power_law_noise)) ./ std(power_law_noise)
    return power_law_noise
end