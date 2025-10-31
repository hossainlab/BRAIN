# Assets Directory

This directory contains images, logos, and other static assets for the BRAIN course website.

## Current Assets

### Logos and Branding
- `logo.svg` - Main BRAIN course logo with stylized brain and text (200x200px)
- `favicon.svg` - Browser favicon with simplified brain icon (32x32px)

Both logos feature:
- Neuroscience-inspired brain illustration
- Purple to pink gradient (#8B5CF6 to #EC4899)
- Neural network connection points
- "BRAIN" text branding

### Logo Design Details

The logo uses a stylized brain design with:
- **Left & Right Hemispheres**: Representing the dual nature of the brain
- **Neural Connections**: Small circles and lines symbolizing synapses
- **Color Gradient**: Purple to pink gradient matching the website theme
- **Typography**: Bold "BRAIN" text with subtitle

## Additional Assets (Optional)

You can add:
- Course promotional images
- Instructor photos
- Diagrams and illustrations
- Banner images
- Module-specific graphics

## Image Guidelines

- **Format**: SVG preferred for logos/icons, PNG or JPEG for photos
- **Size**: Optimize images for web (compress without losing quality)
- **Naming**: Use descriptive, lowercase names with hyphens (e.g., `course-banner.png`)
- **Colors**: Match the website theme (purples: #8B5CF6, pinks: #EC4899)

## Customizing the Logo

To replace the logo with your own design:

1. Create your logo as SVG, PNG, or JPG
2. Save it in this directory as `logo.svg` (or `logo.png`)
3. Update `_quarto.yml` if using different format:
   ```yaml
   navbar:
     logo: assets/your-logo.png
   ```

## Logo Usage

The logo is displayed in:
- Navigation bar (top left)
- Browser tab (favicon)
- Social media previews (when shared)
