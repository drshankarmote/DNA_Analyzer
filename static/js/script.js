// JavaScript to toggle dark mode
function toggleTheme() {
    // Toggle a class that changes the theme (like dark/light mode)
    document.body.classList.toggle('dark-theme');

    // Optionally, you can save the user's theme preference in localStorage
    if (document.body.classList.contains('dark-theme')) {
        localStorage.setItem('theme', 'dark');
    } else {
        localStorage.setItem('theme', 'light');
    }

    function toggleComparisonDetails() {
        var details = document.getElementById('comparisonDetails');
        var button = document.querySelector('button');
        if (details.style.display === "none") {
            details.style.display = "block";
            button.innerHTML = "Hide Detailed Comparison";
        } else {
            details.style.display = "none";
            button.innerHTML = "Show Detailed Comparison";
        }
    }
    
}

// Apply the saved theme preference on page load
window.onload = function () {
    const savedTheme = localStorage.getItem('theme');
    if (savedTheme === 'dark') {
        document.body.classList.add('dark-theme');
    }
}
