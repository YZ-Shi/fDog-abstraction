(function(root) {
    'use strict';

    // A function that applies a flow-based Difference of Gaussians edge-detection filter,
    // based on Kang et al. (2009)
    function fDogFilter(image, value) {

        const kernel = 5;
        let iteration = 2;
        //const sigma_m = 1.0, sigma_c = 0.5;
        const sigma_m = 3.0, sigma_c = 1.0;

        // initialize flow chart
        let result = sobelKernel(image);
        let flow = result[0], gradientMag = result[1];
        flow = rotateFlow(flow, 25.0);
        for (let i = 0; i < iteration; i++) {
            flow = refineFlow(image, flow, gradientMag, kernel);
        }

        image = fDoG(image, flow, sigma_m, sigma_c);
        return image;
    }

    // Construct the flow field
    function flowField(image, value) {
        let mu = 1.5;
        // optional: gaussian blurred

        let newImg = image.createImg(image.width, image.height);

        for (let x = 0; x < image.width; x++) {
            for (let y = 0; y < image.height; y++) {
                let pixel = image.getPixel(x, y);
                let newPx = new Pixel(0, 0, 0);
                let weights = new Array(2 * winR + 1);
                for (let i = 0; i < weights.length; i++) {
                    weights[i] = new Array(2 * winR + 1);
                }
                let wSum = 0;

                for (let i = -winR; i <= winR; i++) {
                    for (let j = -winR; j <= winR; j++) {
                        // border conditions
                        let newX, newY;
                        if (x + i < 0) newX = image.width + x + i;
                        else if (x + i >= image.width) newX = x + i - image.width;
                        else newX = x + i;
                        if (y + j < 0) newY = image.height + y + j;
                        else if (y + j >= image.height) newY = y + j - image.height;
                        else newY = y + j;

                        let neighbor = image.getPixel(newX, newY);
                        let colorDist = Math.sqrt((pixel.data[0] * (l - 1) - neighbor.data[0] * (l - 1)) ** 2 + (pixel.data[1] * (l - 1) - neighbor.data[1] * (l - 1)) ** 2 + (pixel.data[2] * (l - 1) - neighbor.data[2] * (l - 1)) ** 2);
                        let w = Math.exp(-(i ** 2 + j ** 2) / (2 * sigmaS ** 2) - (colorDist ** 2) / (2 * sigmaR ** 2));
                        weights[i + winR][j + winR] = w;
                        wSum += w;
                        }
                }
                for (let i = -winR; i <= winR; i++) {
                    for (let j = -winR; j <= winR; j++) {
                        weights[i + winR][j + winR] /= wSum;

                        let newX, newY;
                        if (x + i < 0) newX = image.width + x + i;
                        else if (x + i >= image.width) newX = x + i - image.width;
                        else newX = x + i;
                        if (y + j < 0) newY = image.height + y + j;
                        else if (y + j >= image.height) newY = y + j - image.height;
                        else newY = y + j;

                        let neighbor = image.getPixel(newX, newY);
                        newPx = newPx.plus(neighbor.multipliedBy(weights[i + winR][j + winR]));
                    }
                }
                newImg.setPixel(x, y, newPx);
            }
        }
        image = newImg;

    }
    function grayscaleFilter(image) {
    for (let x = 0; x < image.width; x++) {
        for (let y = 0; y < image.height; y++) {
            const pixel = image.getPixel(x, y);
            const luminance = 0.2126 * pixel.data[0] + 0.7152 * pixel.data[1] + 0.0722 * pixel.data[2];
            pixel.data[0] = luminance;
            pixel.data[1] = luminance;
            pixel.data[2] = luminance;

            image.setPixel(x, y, pixel);
        }
    }

    return image;
};

    // Compute gradient map with sobel operation
    function gradientMap(image) {
        // discard rgb information
        image = grayscaleFilter(image);
        let imageData = image.getImageData();

        const result = Sobel(imageData);
        console.log(result);
        const sobelData = result;

        const sobelImg = new Image(image.width, image.height, sobelData);
        return sobelImg;
    }

    function sobelKernel(image) {
        var kernelX = [
            [-1,0,1],
            [-2,0,2],
            [-1,0,1]
        ];
    
        var kernelY = [
            [-1,-2,-1],
            [0,0,0],
            [1,2,1]
        ];

        image = Filters.grayscaleFilter(image);

        var grad_x = new Array(image.width), 
        grad_y = new Array(image.width);
        let gradientMag = new Array(image.width),
        flow = new Array(image.width);

        for (let i = 0; i < image.width; i++) {
        grad_x[i] = new Array(image.height).fill(0);
        grad_y[i] = new Array(image.height).fill(0);
        gradientMag[i] = new Array(image.height).fill(0);
        flow[i] = new Array(image.height).fill(0);
        }

        for (let x = 0; x < image.width; x++) {
            for (let y = 0; y < image.height; y++) {
                const pixelX = convolute(image, kernelX, x, y);
                grad_x[x][y] = pixelX.data[0];
                const pixelY = convolute(image, kernelY, x, y);
                grad_y[x][y] = pixelY.data[0];
                let magnitude = Math.sqrt(pixelX.data[0] ** 2 + pixelY.data[0] ** 2);
                gradientMag[x][y] = magnitude;
                let temp = {x: pixelX.data[0]/magnitude, y: pixelY.data[0]/magnitude};
                if (magnitude === 0) flow[x][y] = {x: pixelX.data[0], y: pixelY.data[0]};
                else flow[x][y] = {x: pixelX.data[0]/magnitude, y: pixelY.data[0]/magnitude};
            }
        }
        return [flow, gradientMag];
    }

    function rotateFlow(flow, theta) {
        theta = theta / 180 * pi;
        for(let x = 0; x < flow.length; x++) {
            for (let y = 0; y < flow[0].length; y++) {
                let v = flow[x][y];
                let rx = v.x * Math.cos(theta) - v.y * Math.sin(theta);
                let ry = v.y * Math.cos(theta) + v.x * Math.sin(theta);
                flow[x][y] = {x: rx, y: ry};
            }
        }
        return flow;
    }

    function refineFlow(image, flow, gradientMag, kernel) {
        let refined = new Array(flow.length);
        for (let i = 0; i < flow.length; i++) {
            refined[i] = new Array(flow[0].length).fill(0);
        }
        for (let x = 0; x < image.width; x++) {
            for (let y = 0; y < image.height; y++) {
                refined[x][y] = computeNewVector(image, flow, gradientMag, x, y, kernel, "row");
            }
        }
        flow = refined;
        for (let x = 0; x < image.width; x++) {
            for (let y = 0; y < image.height; y++) {
                refined[x][y] = computeNewVector(image, flow, gradientMag, x, y, kernel, "column");
            }
        }
        flow = refined;
        return flow;
    }

    function computeNewVector(image, flow, gradientMag, x, y, kernel, orientation) {
        let t_x = flow[x][y];
        let t_new = {x: 0, y: 0};
        let wSum = 0;
        if (orientation == "row") {
            for (let r = x - kernel; r <= x + kernel; r++) {
                if (r < 0 || r >= image.width) continue;
                let t_y = flow[r][y];
                let phi = computePhi(t_x, t_y);
                let a = {x: x, y: y}, b = {x: r, y: y};
                let w_s = computeWs(a, b, kernel);
                let w_m = computeWm(gradientMag[x][y], gradientMag[r][y]);
                let w_d = computeWd(t_x, t_y);
                let weight = phi * w_s * w_m * w_d;
                t_new = vPlus(t_new, vMultiply(t_y, weight));
                wSum += weight;
            }
        } else {
            for (let c = y - kernel; c <= y + kernel; c++) {
                if (c < 0 || c >= image.height) continue;
                let t_y = flow[x][c];
                let phi = computePhi(t_x, t_y);
                let a = {x: x, y: y}, b = {x: x, y: c};
                let w_s = computeWs(a, b, kernel);
                let w_m = computeWm(gradientMag[x][y], gradientMag[x][c]);
                let w_d = computeWd(t_x, t_y);
                let weight = phi * w_s * w_m * w_d;
                t_new = vPlus(t_new, vMultiply(t_y, weight));
                wSum += weight;
            }
        }
        let magnitude = Math.sqrt(t_new.x ** 2 + t_new.y ** 2);
        t_new = vDivide(t_new, magnitude);
        return t_new;
    }

    function computePhi(t_x, t_y) {
        if (dot(t_x, t_y) > 0) return 1;
        else return -1;
    }

    function computeWs(a, b, r) {
        if (dFromCenter(b.x, b.y, a) < r) return 1;
        else return 0;
    }

    function computeWm(gradmag_x, gradmag_y) {
        let wm = (1 + Math.tanh(gradmag_y - gradmag_x)) / 2;
        return wm;
    }

    function computeWd(t_x, t_y) {
        return Math.abs(dot(t_x, t_y));
    }

    // sigma_m: determines winR along flow; sigma_c: determines winR along gradient
    function fDoG(image, flow, sigma_m, sigma_c) {
        //const sigma_s = 1.6 * sigma_c;
        const sigma_s = 3 * sigma_c;
        const winR_m = Math.round(sigma_m * 3), winR_c = Math.round(sigma_c * 3), winR_s = Math.round(sigma_s * 3);
        const rho = 0.99; // controls noises
        let h_g = new Array(flow.length), h_e = new Array(flow.length);
        for (let i = 0; i < flow.length; i++) {
            h_g[i] = new Array(flow[0].length).fill(0); // temp value for phase 1
            h_e[i] = new Array(flow[0].length).fill(0); // final value after phase 2
        }

        for (let x = 0; x < image.width; x++) {
            for (let y = 0; y < image.height; y++) {
                let fl = flow[x][y]; // flow vector
                let gr = {x: -fl.y, y: fl.x}; // gradient vector; perpendicular to flow vector

                // Phase 1: 1D DoG filter along the gradient direction
                for (let i = -winR_m; i <= winR_m; i++) {
                    let sampleX1 = Math.round(x + i * fl.x),
                    sampleY1 = Math.round(y + i * fl.y);
                    if (sampleX1 < 0 || sampleX1 >= image.width || sampleY1 < 0 || sampleY1 >= image.height) continue;
                    let g_c = 0, g_s = 0; // sum for DoG
                    let gSum_c = 0, gSum_s = 0; // accumulative weight
                    
                    // compute g_c
                    for (let j = -winR_c; j <= winR_c; j++) {
                        let sampleX11 = Math.round(sampleX1 + j * gr.x),
                        sampleY11 = Math.round(sampleY1 + j * gr.y);
                        if (sampleX11 < 0 || sampleX11 >= image.width || sampleY11 < 0 || sampleY11 >= image.height) continue;
                        let dist = dFromCenter(sampleX11, sampleY11, {x: sampleX1, y: sampleY1});
                        let g = Math.exp(- (dist ** 2) / (2 * sigma_c ** 2)) / Math.sqrt(2 * pi * sigma_c ** 2);
                        let pixel = Filters.samplePixel(image, Math.round(sampleX11), Math.round(sampleY11), "point");
                        g_c += pixel.data[0] * g;
                        gSum_c += g;
                    }
                    // compute g_s
                    for (let j = -winR_s; j <= winR_s; j++) {
                        let sampleX11 = Math.round(sampleX1 + j * gr.x),
                        sampleY11 = Math.round(sampleY1 + j * gr.y);
                        if (sampleX11 < 0 || sampleX11 >= image.width || sampleY11 < 0 || sampleY11 >= image.height) continue;
                        let dist = dFromCenter(sampleX11, sampleY11, {x: sampleX1, y: sampleY1});
                        let g = Math.exp(- (dist ** 2) / (2 * sigma_s ** 2)) / Math.sqrt(2 * pi * sigma_s ** 2);
                        let pixel = Filters.samplePixel(image, Math.round(sampleX11), Math.round(sampleY11), "point");
                        g_s += pixel.data[0] * g;
                        gSum_s += g;
                    }
                    h_g[Math.round(sampleX1)][Math.round(sampleY1)] = g_c / gSum_c - rho * g_s / gSum_s;
                }

                // Phase 2: 1D Gaussian filter along the flow direction
                let g_m = 0, gSum_m = 0;
                for (let i = -winR_m; i <= winR_m; i++) {
                    let sampleX1 = Math.round(x + i * fl.x),
                    sampleY1 = Math.round(y + i * fl.y);
                    if (sampleX1 < 0 || sampleX1 >= image.width || sampleY1 < 0 || sampleY1 >= image.height) continue;
                    let dist = dFromCenter(sampleX1, sampleY1, {x: x, y: y});
                    let g = Math.exp(- (dist ** 2) / (2 * sigma_m ** 2)) / Math.sqrt(2 * pi * sigma_m ** 2);
                    g_m += h_g[sampleX1][sampleY1] * g;
                    gSum_m += g;
                }
                h_e[x][y] = g_m / gSum_m;
            }
        }
        const newImg = image.createImg(image.width, image.height);
        for (let x = 0; x < image.width; x++) {
            for (let y = 0; y < image.height; y++) {
                if (h_e[x][y] < 0 && (1 + Math.tanh(h_e[x][y])) < 1)
                newImg.setPixel(x, y, new Pixel(0, 0, 0));
                else newImg.setPixel(x, y, new Pixel(1, 1, 1));
            }
        }
        return newImg;
    }

})(this);