/**
 * Stitch-inspired Animations for femlabpy Documentation
 * 
 * Uses Intersection Observer to add fade-in-up animations 
 * to cards, sections, and prominent blocks as they scroll into view.
 */

document.addEventListener('DOMContentLoaded', () => {
    // Selectors for elements we want to animate
    const animateSelectors = [
        '.fp-card',
        '.fp-stat',
        '.fp-section-note',
        '.admonition',
        'table.docutils',
        '.bd-article h2',
        '.bd-article h3',
        '.math.notranslate',
        'div.highlight'
    ];

    // Find all matching elements
    const elementsToAnimate = document.querySelectorAll(
        animateSelectors.join(', ')
    );

    // Initial state: hide elements and add transition classes
    elementsToAnimate.forEach((el) => {
        el.classList.add('fp-animate-hidden');
    });

    // Configuration for the Intersection Observer
    const observerOptions = {
        root: null,
        rootMargin: '0px 0px -50px 0px', // Trigger slightly before the element enters the viewport
        threshold: 0.1 // Trigger when 10% of the element is visible
    };

    // The Intersection Observer callback
    const observerCallback = (entries, observer) => {
        entries.forEach((entry) => {
            if (entry.isIntersecting) {
                // Add the visible class to trigger the CSS animation
                entry.target.classList.remove('fp-animate-hidden');
                entry.target.classList.add('fp-animate-visible');
                
                // Stop observing once animated to keep it visible
                observer.unobserve(entry.target);
            }
        });
    };

    // Initialize the observer
    const observer = new IntersectionObserver(observerCallback, observerOptions);

    // Start observing elements
    elementsToAnimate.forEach((el) => observer.observe(el));
});