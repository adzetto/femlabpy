/**
 * Stitch-inspired Animations for femlabpy Documentation
 * 
 * Uses Intersection Observer to add fade-in-up animations 
 * to cards, sections, and prominent blocks as they scroll into view.
 */

document.addEventListener('DOMContentLoaded', () => {
    // Scroll-to-top button
    const topBtn = document.createElement('button');
    topBtn.className = 'fp-scroll-top';
    topBtn.setAttribute('aria-label', 'Scroll to top');
    topBtn.innerHTML = '\u2191';
    document.body.appendChild(topBtn);
    topBtn.addEventListener('click', () => window.scrollTo({ top: 0, behavior: 'smooth' }));
    window.addEventListener('scroll', () => {
        topBtn.classList.toggle('visible', window.scrollY > 400);
    }, { passive: true });

    // Copy-to-clipboard for command rows
    document.querySelectorAll('.fp-cmd-copy').forEach(btn => {
        btn.addEventListener('click', () => {
            const targetId = btn.getAttribute('data-target');
            const codeEl = document.getElementById(targetId);
            if (!codeEl) return;
            const text = codeEl.textContent.trim();
            navigator.clipboard.writeText(text).then(() => {
                btn.classList.add('copied');
                const origHTML = btn.innerHTML;
                btn.innerHTML = '\u2713';
                setTimeout(() => { btn.classList.remove('copied'); btn.innerHTML = origHTML; }, 1500);
            });
        });
    });

    // Selectors for elements we want to animate
    const animateSelectors = [
        '.fp-hero',
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

    // ----------------------------------------------------
    // Hero Naming Animation Logic
    // ----------------------------------------------------
    const animContainer = document.getElementById('fp-naming-anim');
    if (animContainer) {
        const animations = [
            {
                prefix: "k", prefixDesc: "Stiffness Matrix <i>(m for Mass, f for Force)</i>", prefixColor: "#ea4335",
                body: "t3", bodyDesc: "3-Node Triangle <i>(q4 for Quad, h8 for Hex)</i>", bodyColor: "#1a73e8",
                suffix: "e", suffixDesc: "Element-Level Routine", suffixColor: "#34a853",
                argsHtml: `<span class="fp-token" style="color: #fbbc04">X</span>, <span class="fp-token" style="color: #fbbc04">G</span>`,
                argsList: [
                    { name: "X", desc: "Nodal Coordinates", color: "#fbbc04" },
                    { name: "G", desc: "Integration Points & Weights", color: "#fbbc04" }
                ]
            },
            {
                prefix: "m", prefixDesc: "Mass Matrix <i>(k for Stiffness, f for Force)</i>", prefixColor: "#ea4335",
                body: "q4", bodyDesc: "4-Node Quad <i>(t3 for Tri, h8 for Hex)</i>", bodyColor: "#1a73e8",
                suffix: "e", suffixDesc: "Element-Level Routine", suffixColor: "#34a853",
                argsHtml: `<span class="fp-token" style="color: #fbbc04">M</span>, <span class="fp-token" style="color: #fbbc04">T</span>, <span class="fp-token" style="color: #fbbc04">X</span>, <span class="fp-token" style="color: #fbbc04">G</span>`,
                argsList: [
                    { name: "M", desc: "Material Properties Table", color: "#fbbc04" },
                    { name: "T", desc: "Topology Connectivity", color: "#fbbc04" },
                    { name: "X", desc: "Nodal Coordinates", color: "#fbbc04" },
                    { name: "G", desc: "Quadrature Data", color: "#fbbc04" }
                ]
            },
            {
                prefix: "f", prefixDesc: "Force Vector", prefixColor: "#ea4335",
                body: "h8", bodyDesc: "8-Node Hexahedron", bodyColor: "#1a73e8",
                suffix: "e", suffixDesc: "Element-Level Routine", suffixColor: "#34a853",
                argsHtml: `<span class="fp-token" style="color: #fbbc04">F</span>, <span class="fp-token" style="color: #fbbc04">T</span>, <span class="fp-token" style="color: #fbbc04">X</span>, <span class="fp-token" style="color: #fbbc04">G</span>`,
                argsList: [
                    { name: "F", desc: "Face Loads / Body Forces", color: "#fbbc04" },
                    { name: "T", desc: "Topology Connectivity", color: "#fbbc04" },
                    { name: "X", desc: "Nodal Coordinates", color: "#fbbc04" },
                    { name: "G", desc: "Quadrature Data", color: "#fbbc04" }
                ]
            },
            {
                prefix: "q", prefixDesc: "Post-Processing (Internal Forces/Stresses)", prefixColor: "#ea4335",
                body: "t4", bodyDesc: "4-Node Tetrahedron", bodyColor: "#1a73e8",
                suffix: "e", suffixDesc: "Element-Level Routine", suffixColor: "#34a853",
                argsHtml: `<span class="fp-token" style="color: #fbbc04">q</span>, <span class="fp-token" style="color: #fbbc04">T</span>, <span class="fp-token" style="color: #fbbc04">X</span>, <span class="fp-token" style="color: #fbbc04">G</span>, <span class="fp-token" style="color: #fbbc04">u</span>`,
                argsList: [
                    { name: "q", desc: "Material State / History", color: "#fbbc04" },
                    { name: "T, X, G", desc: "Element Geometry & Integ.", color: "#fbbc04" },
                    { name: "u", desc: "Global Displacements", color: "#fbbc04" }
                ]
            }
        ];

        let currentIndex = 0;
        const argPanel = document.getElementById('fp-arg-panel');

        function renderAnimation() {
            const data = animations[currentIndex];

            // Trigger container animation reflow
            animContainer.style.animation = 'none';
            void animContainer.offsetHeight;
            animContainer.style.animation = 'fp-container-cycle 10s forwards';

            // Left side: code window + naming legend
            animContainer.innerHTML = `
                <div class="fp-code-window">
                    <span class="fp-token fp-anim-token" style="--anim-color: ${data.prefixColor}; animation-delay: 0s;">${data.prefix}</span><span class="fp-token fp-anim-token" style="--anim-color: ${data.bodyColor}; animation-delay: 1.5s;">${data.body}</span><span class="fp-token fp-anim-token" style="--anim-color: ${data.suffixColor}; animation-delay: 3s;">${data.suffix}</span><span class="fp-token">(</span><span class="fp-anim-args" style="animation-delay: 4.5s;">${data.argsHtml}</span><span class="fp-token">)</span>
                </div>
                <div class="fp-math-defs">
                    <div class="fp-math-def fp-anim-row" style="animation-delay: 0s;">
                        <span class="fp-math-symbol" style="color: ${data.prefixColor}">${data.prefix}</span>
                        <span class="fp-math-arrow">\u27F6</span>
                        <span class="fp-math-desc">${data.prefixDesc}</span>
                    </div>
                    <div class="fp-math-def fp-anim-row" style="animation-delay: 1.5s;">
                        <span class="fp-math-symbol" style="color: ${data.bodyColor}">${data.body}</span>
                        <span class="fp-math-arrow">\u27F6</span>
                        <span class="fp-math-desc">${data.bodyDesc}</span>
                    </div>
                    <div class="fp-math-def fp-anim-row" style="animation-delay: 3s;">
                        <span class="fp-math-symbol" style="color: ${data.suffixColor}">${data.suffix}</span>
                        <span class="fp-math-arrow">\u27F6</span>
                        <span class="fp-math-desc">${data.suffixDesc}</span>
                    </div>
                </div>
            `;

            // Right side: static arg cards (visible immediately)
            if (argPanel) {
                let cardsHtml = `<div class="fp-arg-header">Arguments</div>`;
                data.argsList.forEach(arg => {
                    cardsHtml += `
                        <div class="fp-arg-card">
                            <span class="fp-arg-badge">${arg.name}</span>
                            <span class="fp-arg-text">${arg.desc}</span>
                        </div>
                    `;
                });
                argPanel.innerHTML = cardsHtml;
            }

            currentIndex = (currentIndex + 1) % animations.length;
        }

        renderAnimation();
        setInterval(renderAnimation, 10000);
    }
});