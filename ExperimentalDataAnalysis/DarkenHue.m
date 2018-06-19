function lighterRGB = DarkenHue(inputRGB, varargin)

    if nargin == 2;
        hueBump = varargin{1};
    else
        hueBump = 0.6;
    end

    hsvHere = rgb2hsv(inputRGB);
    hsvHere(3) = hsvHere(3)*hueBump;
    lighterRGB = hsv2rgb(hsvHere);

end