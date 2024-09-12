import cv2
import numpy as np
import matplotlib.pyplot as plt

# Load the image
image_path = '/content/Image.jpg'  # Replace with the correct path to your image
image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)  # Ensure the image is read as grayscale

# Check if the image was loaded correctly
if image is None:
    print("Error: Image not loaded correctly. Check the file path.")
else:
    # Display the image for reference
    plt.figure(figsize=(10, 8))
    plt.imshow(image, cmap='gray')
    plt.title('Gel Image')
    plt.show()

    # Step 1: Manually define the y-coordinates of the marker bands
    # Example values below, you can adjust these after visually inspecting the image
    y_coords = {
        250: 10,    # 250 kDa at pixel 10 from top
        150: 40,    # 150 kDa at pixel 40 from top
        100: 70,    # 100 kDa at pixel 70 from top
        75: 100,    # 75 kDa at pixel 100 from top
        50: 130,    # 50 kDa at pixel 130 from top
        37: 160,    # 37 kDa at pixel 160 from top
        25: 190,    # 25 kDa at pixel 190 from top
        20: 220,    # 20 kDa at pixel 220 from top
        15: 250     # 15 kDa at pixel 250 from top
    }

    # Step 2: Convert molecular weights to their logarithmic scale
    log_molecular_weights = np.log10(list(y_coords.keys()))
    pixel_positions = np.array(list(y_coords.values()))

    # Step 3: Plot the standard curve (log(MW) vs pixel position)
    plt.figure(figsize=(8, 6))
    plt.scatter(pixel_positions, log_molecular_weights, color='blue', label='Marker Bands')
    plt.title('Standard Curve: log(MW) vs Pixel Position')
    plt.xlabel('Pixel Position')
    plt.ylabel('log(Molecular Weight)')
    plt.grid(True)
    plt.legend()
    plt.show()

    # Step 4: Fit a linear model to get the relationship between pixel position and log(MW)
    coefficients = np.polyfit(pixel_positions, log_molecular_weights, 1)
    print("Fit coefficients (slope, intercept):", coefficients)

    # Step 5: Use the linear model to predict molecular weight of unknown bands
    def predict_molecular_weight(pixel_pos, coefficients):
        log_mw = np.polyval(coefficients, pixel_pos)
        return 10 ** log_mw  # Convert back from log scale to kDa

    # Example: Predict molecular weight of a band at pixel position 180
    unknown_band_pixel_pos = 180
    predicted_mw = predict_molecular_weight(unknown_band_pixel_pos, coefficients)
    print(f"Predicted molecular weight of band at pixel position {unknown_band_pixel_pos}: {predicted_mw:.2f} kDa")
