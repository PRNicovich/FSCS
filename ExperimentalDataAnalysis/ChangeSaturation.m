function lighterRGB = ChangeSaturation(inputRGB, satChange)

    hsvHere = rgb2hsv(inputRGB);
    hsvHere(2) = hsvHere(2)*satChange;
    lighterRGB = hsv2rgb(hsvHere);

end