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

function power_law_noise(len::Integer, NoiseType)
    white_noise = randn(len)
    white_noise_fft = FFTW.rfft(white_noise)
    frequencies = FFTW.rfftfreq(len)
    # Adjust the amplitude of each frequency component according to the power law
    # Avoid division by zero for the first element
    power_law_spectrum = white_noise_fft .* (frequencies .+ 1e-10).^ (NoiseType / 2)
    # Inverse FFT to get noise in time domain
    noises = real(FFTW.irfft(power_law_spectrum, len))
    # Center and normalize the noise
    noises = (noises .- mean(noises)) ./ std(noises)
    return noises
end