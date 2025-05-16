// Apply the saved theme preference on page load
window.onload = function () {
    const toggle = document.getElementById('themeToggle');
    const icon = document.querySelector('.slider .icon');
    const savedTheme = localStorage.getItem('theme');

    if (savedTheme === 'dark') {
        document.body.classList.add('dark-theme');
        toggle.checked = true;
        icon.textContent = '🌙';
    }

    toggle.addEventListener('change', function () {
        document.body.classList.toggle('dark-theme');
        const isDark = document.body.classList.contains('dark-theme');
        icon.textContent = isDark ? '🌙' : '🌞';
        localStorage.setItem('theme', isDark ? 'dark' : 'light');
    });
};


// Optional: Expand/collapse comparison details
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
