#include <stdio.h>  // printf, scanf, fopen, fclose, getchar, putchar
#include <stdlib.h> // standard library header
#include <stdint.h> // intruduces fixed width integer
#include <math.h>   // Used for fabs()  absolute value for float
#include <string.h>
#define M_PI 3.14159265358979323846

typedef struct
{
    uint8_t blue;
    uint8_t green;
    uint8_t red;
} PIXEL;

// p_offset stores the address from where the actual pixel data starts
uint32_t p_offset, width, height, rowSize, rowSize_wo_padding, totalSize, padded_bytes;

// Function for reading an image
uint8_t *Read_Image(FILE *image)
{

    uint32_t bytePerPixel = 3;
    rowSize_wo_padding = width * bytePerPixel;
    rowSize = (width * bytePerPixel + 3) & ~3;
    padded_bytes = rowSize - rowSize_wo_padding;

    totalSize = rowSize * height;

    fseek(image, p_offset, SEEK_SET);
    uint8_t *data = malloc(totalSize);
    fread(data, 1, totalSize, image);

    return data;
}

// function for getting the PIXEL VALUE

PIXEL getPixelValue(uint8_t *imagedata, uint8_t row, uint8_t column)
{
    uint32_t bytePerPixel = 3;

    uint8_t *px = imagedata + (rowSize * row) + (column * bytePerPixel);
    PIXEL *PIXELDATA = (PIXEL *)px;

    return *PIXELDATA;
}

// 3) function for performing quantized value of grayScale image Pixel

uint8_t quantizedVal(PIXEL px, uint8_t level)
{

    // here calculate the original grayScale image PIXEL intensity
    uint8_t gray = (uint8_t)(0.299 * px.red + 0.587 * px.green + 0.114 * px.blue);

    // quantized to reduce the grayScale value
    int step = 255 / (level - 1);
    int newGray = (gray / step) * step;

    return (uint8_t)newGray;
}

// 3) function for writing the grayscale image to new file
uint8_t writeImage(const char *fileName, FILE *srcImage, uint8_t *imageData)
{

    fseek(srcImage, 0, SEEK_SET);
    uint8_t *header = malloc(p_offset);
    fread(header, 1, p_offset, srcImage);

    // open a new file for writing
    FILE *output = fopen(fileName, "wb");
    if (!output)
    {
        printf("Error during creation of new File");
        // Here we have to free the size allocated to header
        free(header);
        return -1;
    }

    // write header to new file
    fwrite(header, 1, p_offset, output);

    // write modified pixel
    fwrite(imageData, 1, totalSize, output);

    fclose(output);
    free(header);

    return 0;
}

// 4) Crop the image

// 4.1) crop function :- return a uint8_t pointer to memory block with image data of crop image

uint8_t *Crop_Image(uint8_t *imagedata, uint32_t startRow, uint32_t startCol, uint32_t cropWidth, uint32_t cropHeight)
{
    uint32_t bytePerPixel = 3;

    // Original row size with padding
    uint32_t srcRowSize = rowSize;

    // new row size with padding
    uint32_t destRowSize = (cropWidth * bytePerPixel + 3) & ~3;
    uint8_t *crop_data = malloc(destRowSize * cropHeight);

    for (uint32_t row = 0; row < cropHeight; row++)
    {
        uint8_t *srcPtr = imagedata + ((startRow + row) * srcRowSize) + (startCol * bytePerPixel);
        uint8_t *destPtr = crop_data + (row * destRowSize);

        // copy actual pixels
        for (uint32_t col = 0; col < cropWidth * bytePerPixel; col++)
        {
            destPtr[col] = srcPtr[col];
        }

        // fill padding bytes with zero          // it is better to fill zero instead than leave a garbage value with it
        for (uint32_t pad = cropWidth * bytePerPixel; pad < destRowSize; pad++)
        {
            destPtr[pad] = 0;
        }
    }

    return crop_data;
}

// 4.2) function to write cropped image(header, data) in file

int writeCroppedImage(const char *fileName, FILE *srcImage, uint8_t *cropData, uint32_t cropWidth, uint32_t cropHeight)
{

    fseek(srcImage, 0, SEEK_SET);
    uint8_t *header = malloc(p_offset);
    fread(header, 1, p_offset, srcImage);

    // writing crop width and height in the crop image header
    *(uint32_t *)(header + 18) = cropWidth;
    *(uint32_t *)(header + 22) = cropHeight;

    // recompute rowSize and total size for cropped image
    uint32_t bytePerPixel = 3;
    uint32_t cropRowSize = (cropWidth * bytePerPixel + 3) & ~3;
    uint32_t cropTotalSize = cropRowSize * cropHeight;

    // Fix file size in header  (offset 2 in BMP header)
    *(uint32_t *)(header + 2) = p_offset + cropTotalSize;

    FILE *output = fopen(fileName, "wb");
    if (!output)
    {
        printf("Error: could not create file \n");
        free(header);
        return -1;
    }

    fwrite(header, 1, p_offset, output);
    fwrite(cropData, 1, cropTotalSize, output);

    fclose(output);
    free(header);

    return 0;
}

// 5) Image enhancement using histogram specification and equalization....

//  5.1) compute histogram
void computeHistogram(uint8_t *imageData, int histogram[256])
{
    for (int i = 0; i < 256; i++)
        histogram[i] = 0;

    for (uint32_t row = 0; row < height; row++)
    {
        for (uint32_t col = 0; col < width; col++)
        {
            PIXEL *px = (PIXEL *)(imageData + (row * rowSize) + (col * 3));
            uint8_t gray = (uint8_t)(0.299 * px->red + 0.587 * px->green + 0.114 * px->blue);
            histogram[gray]++;
        }
    }
}

//  5.2) Histogram Equalization
void histogramEqualization(uint8_t *imagedata)
{
    int hist[256], cdf[256];
    computeHistogram(imagedata, hist);

    int totalPixels = width * height;

    // CDF computation
    cdf[0] = hist[0];
    for (int i = 1; i < 256; i++)
    {
        cdf[i] = cdf[i - 1] + hist[i];
    }

    // Normalize CDF
    uint8_t lookup[256];
    for (int i = 0; i < 256; i++)
    {
        lookup[i] = (uint8_t)(((float)cdf[i] / totalPixels) * 255 + 0.5);
    }

    // Map pixels
    for (uint32_t row = 0; row < height; row++)
    {
        for (uint32_t col = 0; col < width; col++)
        {
            PIXEL *px = (PIXEL *)(imagedata + (row * rowSize) + (col * 3));
            uint8_t gray = (uint8_t)(0.299 * px->red + 0.587 * px->green + 0.114 * px->blue);
            uint8_t newGray = lookup[gray];
            px->red = px->green = px->blue = newGray;
        }
    }
}

//  5.3) Histogram Specification
void histogramSpecification(uint8_t *imagedata, int targetHist[256])
{
    int hist[256], cdf[256], targetCdf[256];
    computeHistogram(imagedata, hist);

    int totalPixels = width * height;

    // compute input CDF
    cdf[0] = hist[0];
    for (int i = 1; i < 256; i++)
        cdf[i] = cdf[i - 1] + hist[i];

    // Normalize input CDF
    float cdfNorm[256];
    for (int i = 0; i < 256; i++)
    {
        cdfNorm[i] = (float)cdf[i] / totalPixels;
    }

    // Target CDF
    targetCdf[0] = targetHist[0];
    for (int i = 1; i < 256; i++)
    {
        targetCdf[i] = targetCdf[i - 1] + targetHist[i];
    }

    float targetCdfNorm[256];
    int targetTotal = targetCdf[255];
    for (int i = 0; i < 256; i++)
    {
        targetCdfNorm[i] = (float)targetCdf[i] / targetTotal;
    }

    // Build mapping table
    uint8_t mapping[256];
    for (int i = 0; i < 256; i++)
    {
        float diff = 1.0;
        int best = 0;
        //      If my input gray level has cumulative probability cdfNorm[i], which target gray level j has the closest cumulative probability?

        //      That’s how we "remap" the intensities so that the shape of the input histogram becomes as close as possible to the target histogram.
        for (int j = 0; j < 256; j++)
        {
            float d = fabs(cdfNorm[i] - targetCdfNorm[j]);
            if (d < diff)
            {
                diff = d;
                best = j;
            }
        }
        mapping[i] = best;
    }

    // Map Pixels
    for (uint32_t row = 0; row < height; row++)
    {
        for (uint32_t col = 0; col < width; col++)
        {
            PIXEL *px = (PIXEL *)(imagedata + (row * rowSize) + (col * 3));
            uint8_t gray = (uint8_t)(0.299 * px->red + 0.587 * px->green + 0.114 * px->blue);
            uint8_t newGray = mapping[gray];
            px->red = px->green = px->blue = newGray;
        }
    }
}

// 6) Rotate image 90degrees or 180degrees

// 6.1) Rotate image 90 degrees clockwise
uint8_t *rotate90(uint8_t *imagedata)
{
    int newWidth = height;
    int newHeight = width;
    int newRowSize = (newWidth * 3 + 3) & (~3); // row padding to multiple of 4
    uint8_t *rotated = (uint8_t *)malloc(newRowSize * newHeight);

    for (uint32_t row = 0; row < height; row++)
    {
        for (uint32_t col = 0; col < width; col++)
        {
            PIXEL *src = (PIXEL *)(imagedata + row * rowSize + col * 3);
            PIXEL *dst = (PIXEL *)(rotated + col * newRowSize + (newWidth - row - 1) * 3);
            *dst = *src;
        }
    }

    // width = newWidth;
    // height = newHeight;
    // rowSize = newRowSize;
    return rotated;
}

// 6.2) Rotate image 180 degrees clockwise
uint8_t *rotate180(uint8_t *imagedata, int width, int height, int rowSize)
{
    uint8_t *rotated = (uint8_t *)malloc(rowSize * height);

    for (uint32_t row = 0; row < height; row++)
    {
        for (uint32_t col = 0; col < width; col++)
        {
            PIXEL *src = (PIXEL *)(imagedata + row * rowSize + col * 3);
            PIXEL *dst = (PIXEL *)(rotated + (height - row - 1) * rowSize + (width - col - 1) * 3);
            *dst = *src;
        }
    }
    return rotated;
}

// 6.3) Write Rotated Image     // function is used for both writing 90deg as well as 180deg
int writeRotatedImage(const char *fileName, FILE *srcImage, uint8_t *rotatedData, uint32_t newWidth, uint32_t newHeight)
{
    // Read original header
    fseek(srcImage, 0, SEEK_SET);
    uint8_t *header = malloc(p_offset);
    fread(header, 1, p_offset, srcImage);

    // Update width and height in the header (offset 18 = width, 22 = height)
    *(uint32_t *)(header + 18) = newWidth;
    *(uint32_t *)(header + 22) = newHeight;

    // Recompute rowSize and total size for rotated image
    uint32_t bytePerPixel = 3;
    uint32_t rotatedRowSize = (newWidth * bytePerPixel + 3) & ~3;
    uint32_t rotatedTotalSize = rotatedRowSize * newHeight;

    // Fix file size in header (offset 2 = file size)
    *(uint32_t *)(header + 2) = p_offset + rotatedTotalSize;

    // Open output file
    FILE *output = fopen(fileName, "wb");
    if (!output)
    {
        printf("Error: could not create rotated file\n");
        free(header);
        return -1;
    }

    // Write header + pixel data
    fwrite(header, 1, p_offset, output);
    fwrite(rotatedData, 1, rotatedTotalSize, output);

    fclose(output);
    free(header);

    return 0;
}

// 6.4) Rotate an image by an arbitrrary angle theta

uint8_t *rotateTheta(uint8_t *imagedata, double angle, int *outWidth, int *outHeight)
{
    double theta = angle * M_PI / 180.0; // convert to radians
    double cosA = cos(theta);
    double sinA = sin(theta);

    // New dimensions after rotation (bounding box)
    int newWidth = (int)(fabs(width * cosA) + fabs(height * sinA));
    int newHeight = (int)(fabs(width * sinA) + fabs(height * cosA));

    // Row size with padding
    int newRowSize = (newWidth * 3 + 3) & ~3;
    uint8_t *rotated = (uint8_t *)calloc(newRowSize * newHeight, 1);

    double cx = width / 2.0;
    double cy = height / 2.0;
    double newCx = newWidth / 2.0;
    double newCy = newHeight / 2.0;

    for (int yNew = 0; yNew < newHeight; yNew++)
    {
        for (int xNew = 0; xNew < newWidth; xNew++)
        {
            // Map destination pixel (xNew,yNew) → source pixel (x,y)
            double xRel = xNew - newCx;
            double yRel = yNew - newCy;

            double x = xRel * cosA + yRel * sinA + cx;
            double y = -xRel * sinA + yRel * cosA + cy;

            int xi = (int)round(x);
            int yi = (int)round(y);

            if (xi >= 0 && xi < (int)width && yi >= 0 && yi < (int)height)
            {
                PIXEL *src = (PIXEL *)(imagedata + yi * rowSize + xi * 3);
                PIXEL *dst = (PIXEL *)(rotated + yNew * newRowSize + xNew * 3);
                *dst = *src;
            }
            // else leave black (calloc already initialized to 0)
        }
    }

    *outWidth = newWidth;
    *outHeight = newHeight;
    return rotated;
}

// 7) Implement contrast enhancement and sharpening of an image.

// 7.1) Implemented contrast enhancement

void contrastEnhancement(uint8_t *imagedata)
{
    uint8_t min = 255, max = 0;

    // Find min and max intensity
    for (uint32_t row = 0; row < height; row++)
    {
        for (uint32_t col = 0; col < width; col++)
        {
            PIXEL *px = (PIXEL *)(imagedata + row * rowSize + col * 3);
            uint8_t gray = (uint8_t)(0.299 * px->red + 0.587 * px->green + 0.114 * px->blue);
            if (gray < min)
                min = gray;
            if (gray > max)
                max = gray;
        }
    }

    // Apply linear contrast stretching
    for (uint32_t row = 0; row < height; row++)
    {
        for (uint32_t col = 0; col < width; col++)
        {
            PIXEL *px = (PIXEL *)(imagedata + row * rowSize + col * 3);

            px->red = (uint8_t)((px->red - min) * 255.0 / (max - min));
            px->green = (uint8_t)((px->green - min) * 255.0 / (max - min));
            px->blue = (uint8_t)((px->blue - min) * 255.0 / (max - min));
        }
    }
}

// 7.2) Implemented sharpening of image

//    Sharpening using 3x3 kernel

void sharpenImage(uint8_t *imagedata)
{
    // Copy original data to avoid overwrite during convolution
    uint8_t *copy = (uint8_t *)malloc(rowSize * height);
    memcpy(copy, imagedata, rowSize * height);

    int kernel[3][3] = {
        {0, -1, 0},
        {-1, 5, -1},
        {0, -1, 0}};

    for (uint32_t row = 1; row < height - 1; row++)
    {
        for (uint32_t col = 1; col < width - 1; col++)
        {
            int r = 0, g = 0, b = 0;

            // Apply convolution
            for (int kr = -1; kr <= 1; kr++)
            {
                for (int kc = -1; kc <= 1; kc++)
                {
                    PIXEL *p = (PIXEL *)(copy + (row + kr) * rowSize + (col + kc) * 3);
                    int kVal = kernel[kr + 1][kc + 1];
                    r += p->red * kVal;
                    g += p->green * kVal;
                    b += p->blue * kVal;
                }
            }

            // Clamp values to 0–255
            if (r < 0)
                r = 0;
            if (r > 255)
                r = 255;
            if (g < 0)
                g = 0;
            if (g > 255)
                g = 255;
            if (b < 0)
                b = 0;
            if (b > 255)
                b = 255;

            PIXEL *dst = (PIXEL *)(imagedata + row * rowSize + col * 3);
            dst->red = (uint8_t)r;
            dst->green = (uint8_t)g;
            dst->blue = (uint8_t)b;
        }
    }

    free(copy);
}

// 8) Bluring the image function 3*3 averaging filter

// 3*3 average filter implementation
void blurImage(uint8_t *imagedata)
{
    // Make a copy of original so we don’t overwrite while processing
    uint8_t *copy = (uint8_t *)malloc(rowSize * height);
    memcpy(copy, imagedata, rowSize * height);

    int kernel[3][3] = {
        {1, 1, 1},
        {1, 1, 1},
        {1, 1, 1}};
    int kernelSum = 9; // sum of kernel values

    for (uint32_t row = 1; row < height - 1; row++)
    {
        for (uint32_t col = 1; col < width - 1; col++)
        {
            int r = 0, g = 0, b = 0;

            // Apply convolution
            for (int kr = -1; kr <= 1; kr++)
            {
                for (int kc = -1; kc <= 1; kc++)
                {
                    PIXEL *p = (PIXEL *)(copy + (row + kr) * rowSize + (col + kc) * 3);
                    r += p->red * kernel[kr + 1][kc + 1];
                    g += p->green * kernel[kr + 1][kc + 1];
                    b += p->blue * kernel[kr + 1][kc + 1];
                }
            }

            r /= kernelSum;
            g /= kernelSum;
            b /= kernelSum;

            // Assign blurred pixel
            PIXEL *dst = (PIXEL *)(imagedata + row * rowSize + col * 3);
            dst->red = (uint8_t)r;
            dst->green = (uint8_t)g;
            dst->blue = (uint8_t)b;
        }
    }

    free(copy);
}

// 5*5 average filter
uint8_t *blur5x5(uint8_t *imagedata)
{
    uint8_t *blurred = (uint8_t *)malloc(rowSize * height);

    int kernelSize = 5;
    int k = kernelSize / 2;                   // = 2 for 5x5
    int kernelArea = kernelSize * kernelSize; // = 25

    for (uint32_t row = 0; row < height; row++)
    {
        for (uint32_t col = 0; col < width; col++)
        {
            float sumR = 0, sumG = 0, sumB = 0;
            int count = 0;

            // Loop over 5x5 neighborhood
            for (int dy = -k; dy <= k; dy++)
            {
                for (int dx = -k; dx <= k; dx++)
                {
                    int nrow = row + dy;
                    int ncol = col + dx;

                    // Check bounds
                    if (nrow >= 0 && nrow < height && ncol >= 0 && ncol < width)
                    {
                        PIXEL *neighbor = (PIXEL *)(imagedata + nrow * rowSize + ncol * 3);
                        sumR += neighbor->red;
                        sumG += neighbor->green;
                        sumB += neighbor->blue;
                        count++;
                    }
                }
            }

            // Take average (use count in case near borders < 25 pixels available)
            PIXEL *dst = (PIXEL *)(blurred + row * rowSize + col * 3);
            dst->red = (uint8_t)(sumR / count);
            dst->green = (uint8_t)(sumG / count);
            dst->blue = (uint8_t)(sumB / count);
        }
    }

    return blurred;
}

// 9) Generate a negative of the image
uint8_t *negativeImage(uint8_t *imagedata)
{
    // Allocate memory for new image (same size as original)
    uint8_t *negative = (uint8_t *)malloc(rowSize * height);

    for (uint32_t row = 0; row < height; row++)
    {
        for (uint32_t col = 0; col < width; col++)
        {
            // Source pixel
            PIXEL *src = (PIXEL *)(imagedata + row * rowSize + col * 3);
            // Destination pixel
            PIXEL *dst = (PIXEL *)(negative + row * rowSize + col * 3);

            // Invert each channel (R, G, B)
            dst->red = 255 - src->red;
            dst->green = 255 - src->green;
            dst->blue = 255 - src->blue;
        }
    }

    return negative;
}

int main()
{

    FILE *image = fopen("image.bmp", "rb"); // After fopen() function some system resources such as buffers are allocated for processing the input file.

    // check whether image is succesfuly loaded or not...
    if (image == NULL)
    {
        printf("Error occured during opening of the file");
        fclose(image); // After the operation on file is over, fclose() deallocate the resources given to the FILE object.
        return 1;
    }

    // mention the position in variabe name

    fseek(image, 10, SEEK_SET);
    fread(&p_offset, sizeof(uint32_t), 1, image);
    fseek(image, 18, SEEK_SET);
    fread(&width, sizeof(uint32_t), 1, image);
    fread(&height, sizeof(uint32_t), 1, image);

    int choice;
    do
    {
        printf("\n===== MENU =====\n");
        printf("1. Get Height & Width\n");
        printf("2. Get Pixel Value\n");
        printf("3. Convert to Quantized GrayScale Image\n");
        printf("4. Crop Image\n");
        printf("5. Histogram Equalization and Histogram Specification\n");
        printf("6. Rotate an image at 90 or 180 or any given angle\n");
        printf("7. Implement contrast enhancement and sharpening of an image. \n");
        printf("8. Bluring the image \n");
        printf("9. generate a negative of an image \n");
        printf("10. Exit\n");
        printf("Enter choice: ");
        scanf("%d", &choice);

        switch (choice)
        {
        case 1:
            printf("Width = %d, Height = %d\n", width, height);
            break;

        case 2:
        {

            // 2) Read pixel value of image
            uint8_t *imagedata = Read_Image(image);

            int r, c;
            printf("Enter row and column: ");
            scanf("%d %d", &r, &c);

            // Get the PIXEL VALUE
            PIXEL px = getPixelValue(imagedata, r, c);
            printf("Pixel Value at (%d,%d) -> R=%d G=%d B=%d\n", r, c, px.red, px.green, px.blue);
            free(imagedata);
            break;
        }

        case 3:
        {

            // 3) reduce the number of gray level of an image (Gray Level quantization)

            // Gray Scale Quantization :- It is the process of reducing the different shades of gray of image into finite distinct set of gray levels..
            // why Gray Scale Quantization is performed ?
            // These is done to reduce the image data, simple storage and processing, and facilitate digital representation by assigning each pixel
            // one of the limited intensity values, often decribed as k distinct level, such as 256 level for 8-bit gray scale image

            // gray level :- basically means intensity of pixel in grayscale image
            // gray scale image does not have colors like red, green, blue instead it contains shades of gray i.e. black and white
            // Each pixel has a numerical value representing how bright/dark the pixel is. These value is called as gray level

            // Gray Level
            // 0 :- pure black   255 :- pure white

            // Function to convert the pixel value to Quantized gray scale value

            // loop for storing the Quantized gray intensity value in the PIXEL datavariables
            uint8_t *imagedata = Read_Image(image);
            for (uint32_t row = 0; row < height; row++)
            {
                for (uint32_t col = 0; col < width; col++)
                {
                    PIXEL *px = (PIXEL *)(imagedata + (row * rowSize) + (col * 3));
                    uint8_t newGray = quantizedVal(*px, 4);
                    px->red = newGray;
                    px->green = newGray;
                    px->blue = newGray;
                }
            }
            if (writeImage("QuantizedGrayScaleImage.bmp", image, imagedata) == 0)
            {
                printf("Quantized GrayScale Image saved!\n");
            }
            free(imagedata);
            break;
        }

        case 4:
        {

            // 4) Crop Image implementation

            uint8_t *imagedata = Read_Image(image);

            // Lets take some points imformation regarding the start of xcordinate, ycordinate of the crop image and cropWidth, cropHeight
            uint32_t xStart, yStart, crop_width, crop_height;
            printf("Enter xStart, yStart, cropWidth, cropHeight: ");
            scanf("%u %u %u %u", &xStart, &yStart, &crop_width, &crop_height);

            if (xStart + crop_width > width || yStart + crop_height > height)
            {
                printf("Crop area is outside the image bound!\n");
                free(imagedata);
                break;
            }

            uint8_t *cropped = Crop_Image(imagedata, yStart, xStart, crop_width, crop_height);

            if (writeCroppedImage("CroppedImage.bmp", image, cropped, crop_width, crop_height) == 0)
            {
                printf("Cropped Image saved as CroppedImage.bmp\n");
            }

            free(imagedata);
            free(cropped);
            break;
        }

        case 5:
        {
            // Image enhancement using histogram specification and equalization.

            uint8_t *imagedata = Read_Image(image);

            // Equalization
            histogramEqualization(imagedata);
            writeImage("HistEqualized.bmp", image, imagedata);

            free(imagedata);

            // Important: reset file pointer before reading original image again
            fseek(image, 0, SEEK_SET);
            imagedata = Read_Image(image);

            // Specification (example: match to dark-biased histogram)
            int targetHist[256];
            for (int i = 0; i < 256; i++)
            {
                targetHist[i] = (i < 64) ? 1 : 10; // strong bright bias
            }
            histogramSpecification(imagedata, targetHist);
            writeImage("HistSpecified.bmp", image, imagedata);

            free(imagedata);
            break;
        }

        case 6:
        {
            uint8_t *imagedata = Read_Image(image);

            int choiceRot;
            printf("1. Rotate 90\n2. Rotate 180\n3. Rotate by custom angle\n");
            scanf("%d", &choiceRot);

            if (choiceRot == 1)
            {
                uint8_t *rotated = rotate90(imagedata);
                writeRotatedImage("Rotated90.bmp", image, rotated, height, width);
                free(rotated);
            }
            else if (choiceRot == 2)
            {
                uint8_t *rotated = rotate180(imagedata, width, height, rowSize);
                writeRotatedImage("Rotated180.bmp", image, rotated, width, height);
                free(rotated);
            }
            else if (choiceRot == 3)
            {
                double angle;
                printf("Enter angle in degrees: ");
                scanf("%lf", &angle);

                int newW, newH;
                uint8_t *rotated = rotateTheta(imagedata, angle, &newW, &newH);
                writeRotatedImage("RotatedTheta.bmp", image, rotated, newW, newH);
                free(rotated);
            }

            free(imagedata);
            break;
        }

        case 7:
        {
            // Implement contrast enhancement and sharpening of an image.

            uint8_t *imagedata = Read_Image(image);

            int subChoice;
            printf("1. Contrast Enhancement\n2. Sharpening\n");
            scanf("%d", &subChoice);

            // 7.1) Contrast Enhancement (Linear Stretching)

            // The idea:
            // Find the minimum and maximum intensity in the image
            // Stretch the range [min, max] → [0, 255]
            // This improves global contrast

            if (subChoice == 1)
            {
                contrastEnhancement(imagedata);
                writeImage("ContrastEnhanced.bmp", image, imagedata);
                printf("Contrast Enhanced Image saved!\n");
            }
            else if (subChoice == 2)
            {
                sharpenImage(imagedata);
                writeImage("Sharpened.bmp", image, imagedata);
                printf("Sharpened Image saved!\n");
            }

            free(imagedata);
            break;
        }

            // contrastEnhancement(imagedata);
            // writeRotatedImage("ContrastEnhanced.bmp", image, imagedata, width, height);
            // free(imagedata);
            // break;

        case 8:
        {
            // Bluring of image using the 3*3 averaging filter(Box Blur)
            // Q) why to blur a image?
            // ans) Images sometimes have sharp noise (random dots, rough texture).We want to smooth or blur the image to reduce that noise.

            // An image is just a big grid of pixels (numbers for brightness or RGB).
            // Instead of taking one pixel’s value as-is, we look at its neighborhood.
            // Then, we replace the pixel with a weighted average of its neighbors.

            // This means:
            // For each pixel, look at itself + its 8 neighbors. Add them all together. Divide by 9 (so it’s the average).

            // 5*5 or 7*7 blur more strongly...

            // Read original image
            uint8_t *imagedata = Read_Image(image);

            // Apply 3x3 blur
            blurImage(imagedata);
            writeImage("Blur3x3.bmp", image, imagedata);

            // Reset file pointer and read again for 5x5 blur
            fseek(image, 0, SEEK_SET);
            uint8_t *imagedata2 = Read_Image(image);

            // Apply 5x5 blur
            uint8_t *blurred5 = blur5x5(imagedata2);
            writeImage("Blur5x5.bmp", image, blurred5);

            free(imagedata);
            free(imagedata2);
            free(blurred5);

            printf("Blurred images saved as Blur3x3.bmp and Blur5x5.bmp\n");
            break;
        }

        case 9:
        {
            // Generate negative of an image
            uint8_t *imagedata = Read_Image(image);

            uint8_t *neg = negativeImage(imagedata);
            writeImage("Negative.bmp", image, neg);

            free(imagedata);
            free(neg);

            printf("Negative image saved as Negative.bmp\n");
            break;
        }

        case 10:
        {
            printf("Exiting...\n");
            break;
        }
        default:
            printf("Invalid choice!\n");
        }
    } while (choice != 10);

    fclose(image);

    return 0;
}
