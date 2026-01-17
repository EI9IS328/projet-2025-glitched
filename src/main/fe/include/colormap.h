//************************************************************************
//   proxy application v.0.0.1
//
//  colormap.h: colormap utilities for visualization
//
//************************************************************************

#ifndef COLORMAP_H_
#define COLORMAP_H_

#include <cstdint>
#include <string>

/**
 * @brief RGB color structure
 */
struct RGB {
  uint8_t r;
  uint8_t g;
  uint8_t b;
};

/**
 * @brief Colormap types
 */
enum ColormapType {
  COLORMAP_GRAYSCALE,  // Use PGM format
  COLORMAP_VIRIDIS,    // Use PPM format
  COLORMAP_JET         // Use PPM format
};

/**
 * @brief Convert a string to a ColormapType
 * @param name The colormap name (case-insensitive)
 * @return The corresponding ColormapType
 */
ColormapType colormapFromString(const std::string& name);

/**
 * @brief Map a normalized value [0-255] to an RGB color using the specified colormap
 * @param value Input value in range [0, 255]
 * @param type Colormap type to use
 * @return RGB color
 */
RGB applyColormap(uint8_t value, ColormapType type);

/**
 * @brief Apply Viridis colormap
 * @param value Input value in range [0, 255]
 * @return RGB color
 */
RGB viridis(uint8_t value);

/**
 * @brief Apply Jet colormap
 * @param value Input value in range [0, 255]
 * @return RGB color
 */
RGB jet(uint8_t value);

#endif /* COLORMAP_H_ */