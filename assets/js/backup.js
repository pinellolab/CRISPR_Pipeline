// Function to open a specific tab and display its content
function openTab(evt, modalityName) {
    var i, tabcontent, tablinks;
    
    // Hide all tab content
    tabcontent = document.getElementsByClassName("tabcontent");
    for (i = 0; i < tabcontent.length; i++) {
        tabcontent[i].style.display = "none";
    }
    
    // Remove 'active' class from all tab links
    tablinks = document.getElementsByClassName("tablinks");
    for (i = 0; i < tablinks.length; i++) {
        tablinks[i].className = tablinks[i].className.replace(" active", "");
        tablinks[i].setAttribute("aria-selected", "false");
    }
    
    // Show the selected tab's content and add 'active' class to the clicked tab
    document.getElementById(modalityName).style.display = "block";
    evt.currentTarget.className += " active";
    evt.currentTarget.setAttribute("aria-selected", "true");
    
    // Display the relevant cards for the selected tab
    displayCards(modalityName);
}

// Function to display a specific set of cards based on the selected tab and page
function displayCards(modalityName, page = 1, cardsPerPage = 9) {
    const tabContent = document.getElementById(modalityName);
    const cards = tabContent.getElementsByClassName('card');
    const totalPages = Math.ceil(cards.length / cardsPerPage);

    // Loop through cards and display only those on the current page
    for (let i = 0; i < cards.length; i++) {
        if (i >= (page - 1) * cardsPerPage && i < page * cardsPerPage) {
            cards[i].style.display = 'block';
        } else {
            cards[i].style.display = 'none';
        }
    }
    
    // Update pagination controls based on the total pages
    updatePaginationControls(modalityName, page, totalPages);
}

// Function to update pagination buttons based on the current page and total pages
function updatePaginationControls(modalityName, currentPage, totalPages) {
    const paginationContainer = document.getElementById(`pagination-${modalityName}`);
    paginationContainer.innerHTML = '';

    for (let i = 1; i <= totalPages; i++) {
        const button = document.createElement('button');
        button.textContent = i;
        button.onclick = () => displayCards(modalityName, i);
        if (i === currentPage) {
            button.classList.add('active');
        }
        paginationContainer.appendChild(button);
    }
}

// Function to filter and sort cards based on user input and selected criteria
function filterAndSortCards() {
    const filterValue = document.getElementById('filter-input').value.toLowerCase();
    const sortValue = document.getElementById('sort-select').value;
    const cards = Array.from(document.querySelectorAll('.card'));

    // Filter cards based on the input value
    cards.forEach(card => {
        const title = card.querySelector('h3').textContent.toLowerCase();
        card.style.display = title.includes(filterValue) ? 'block' : 'none';
    });

    // Sort cards based on the selected criteria
    cards.sort((a, b) => {
        const aTitle = a.querySelector('h3').textContent;
        const bTitle = b.querySelector('h3').textContent;
        const aValue = parseInt(a.querySelector('.highlight').textContent);
        const bValue = parseInt(b.querySelector('.highlight').textContent);

        switch (sortValue) {
            case 'name-asc': return aTitle.localeCompare(bTitle);
            case 'name-desc': return bTitle.localeCompare(aTitle);
            case 'value-asc': return aValue - bValue;
            case 'value-desc': return bValue - aValue;
        }
    });

    // Reattach sorted cards to the container
    const container = document.querySelector('.flex-container');
    cards.forEach(card => container.appendChild(card));
}

// Event listeners for filtering, sorting, and showing images
document.getElementById('filter-input').addEventListener('input', filterAndSortCards);
document.getElementById('sort-select').addEventListener('change', filterAndSortCards);
document.getElementById('show-images').addEventListener('change', function() {
    const images = document.querySelectorAll('.card img');
    images.forEach(img => img.style.display = this.checked ? 'block' : 'none');
});

// Event listener for adjusting the number of cards per row
document.getElementById('cards-per-row').addEventListener('change', function() {
    const cards = document.querySelectorAll('.card');
    const percentage = 100 / parseInt(this.value);
    cards.forEach(card => card.style.width = `calc(${percentage}% - 20px)`);
});

// Function to toggle the display of elements like tables or images

function toggleElement(button, elementType) {
    const card = button.closest('.card');
    let container;
    
    if (elementType === 'table') {
        container = card.querySelector('.table-container');
        const searchContainer = card.querySelector('.table-search-container');
        
        if (container) {
            if (container.style.display === 'none' || container.style.display === '') {
                container.style.display = 'block';
                button.classList.add('active');

                // Show search bar if it exists
                if (searchContainer) {
                    searchContainer.style.display = 'block';
                    const searchInput = searchContainer.querySelector('.table-search');
                    if (searchInput) {
                        // Add event listener for search input
                        searchInput.addEventListener('input', function() {
                            filterTable(container, this.value);
                        });
                    }
                }

                // Initialize pagination when showing the table
                if (!container.dataset.paginationInitialized) {
                    initializeTablePagination(container);
                    container.dataset.paginationInitialized = 'true';
                }
            } else {
                container.style.display = 'none';
                button.classList.remove('active');
                
                // Hide search bar if it exists
                if (searchContainer) {
                    searchContainer.style.display = 'none';
                    // Clear search input and reset table view
                    const searchInput = searchContainer.querySelector('.table-search');
                    if (searchInput) {
                        searchInput.value = '';
                        resetTableSearch(container);
                    }
                }
            }
        }
    } else if (elementType === 'figure') {
        const imgSrc = button.getAttribute('data-imgsrc');
        const modal = document.getElementById('figureModal');
        const modalImg = document.getElementById('modalImage');
        modal.style.display = 'block';
        modalImg.src = imgSrc;
    }

    // Remove active class from all other buttons
    card.querySelectorAll('.toggle-button').forEach(btn => {
        if (btn !== button) {
            btn.classList.remove('active');
        }
    });
}

// Function to initialize search bar for a table
const ROWS_PER_PAGE = 10;

function filterTable(container, searchTerm) {
    const table = container.querySelector('table');
    const rows = table.querySelectorAll('tbody tr');
    let visibleRowCount = 0;

    const searchRegex = new RegExp(searchTerm.replace(/[-\/\\^$*+?.()|[\]{}]/g, '\\$&'), 'i');

    rows.forEach((row, index) => {
        const sgRNA = row.querySelector('td:first-child').textContent;
        if (searchRegex.test(sgRNA)) {
            row.style.display = '';
            row.classList.remove('filtered-out');
            visibleRowCount++;
        } else {
            row.style.display = 'none';
            row.classList.add('filtered-out');
        }
    });

    updatePagination(container, visibleRowCount);
    showPage(container, 1);
}


function updatePagination(container, visibleRowCount) {
    const pageCount = Math.ceil(visibleRowCount / ROWS_PER_PAGE);
    const paginationContainer = container.querySelector('.table-pagination');
    paginationContainer.innerHTML = '';

    if (pageCount > 1) {
        // Previous button
        const prevBtn = document.createElement('button');
        prevBtn.innerText = '<';
        prevBtn.addEventListener('click', () => {
            const currentPage = parseInt(paginationContainer.querySelector('button.active').innerText);
            if (currentPage > 1) {
                showPage(container, currentPage - 1);
            }
        });
        paginationContainer.appendChild(prevBtn);

        // Determine the range of visible pages
        const MAX_VISIBLE_PAGES = 3;
        let currentPage = parseInt(paginationContainer.querySelector('button.active')?.innerText || '1');
        let startPage = Math.max(1, Math.min(currentPage - 2, pageCount - MAX_VISIBLE_PAGES + 1));
        let endPage = Math.min(pageCount, startPage + MAX_VISIBLE_PAGES - 1);

        // First page and ellipsis before visible range if necessary
        if (startPage > 1) {
            const firstPageBtn = document.createElement('button');
            firstPageBtn.textContent = '1';
            firstPageBtn.addEventListener('click', () => showPage(container, 1));
            paginationContainer.appendChild(firstPageBtn);

            const ellipsisStart = document.createElement('span');
            ellipsisStart.innerText = '...';
            paginationContainer.appendChild(ellipsisStart);
        }

        // Visible page buttons
        for (let i = startPage; i <= endPage; i++) {
            const button = document.createElement('button');
            button.textContent = i;
            button.setAttribute('data-page', i);
            button.addEventListener('click', () => showPage(container, i));
            if (i === currentPage) {
                button.classList.add('active');
            }
            paginationContainer.appendChild(button);
        }

        // Ellipsis after visible range and last page if necessary
        if (endPage < pageCount) {
            const ellipsisEnd = document.createElement('span');
            ellipsisEnd.innerText = '...';
            paginationContainer.appendChild(ellipsisEnd);

            const lastPageBtn = document.createElement('button');
            lastPageBtn.textContent = pageCount;
            lastPageBtn.addEventListener('click', () => showPage(container, pageCount));
            paginationContainer.appendChild(lastPageBtn);
        }

        // Next button
        const nextBtn = document.createElement('button');
        nextBtn.innerText = '>';
        nextBtn.addEventListener('click', () => {
            const currentPage = parseInt(paginationContainer.querySelector('button.active').innerText);
            if (currentPage < pageCount) {
                showPage(container, currentPage + 1);
            }
        });
        paginationContainer.appendChild(nextBtn);

        // Direct page entry
        const pageInput = document.createElement('input');
        pageInput.type = 'number';
        pageInput.min = 1;
        pageInput.max = pageCount;
        pageInput.placeholder = 'Go ';
        pageInput.classList.add('pagination-input');
        pageInput.addEventListener('keypress', function(e) {
            if (e.key === 'Enter') {
                const pageNum = parseInt(pageInput.value);
                if (!isNaN(pageNum) && pageNum >= 1 && pageNum <= pageCount) {
                    showPage(container, pageNum);
                }
            }
        });
        paginationContainer.appendChild(pageInput);
    }
}

function showPage(container, pageNum) {
    const table = container.querySelector('table');
    const rows = table.querySelectorAll('tbody tr:not(.filtered-out)');
    const startIndex = (pageNum - 1) * ROWS_PER_PAGE;
    const endIndex = startIndex + ROWS_PER_PAGE;

    rows.forEach((row, index) => {
        if (index >= startIndex && index < endIndex) {
            row.style.display = '';
        } else {
            row.style.display = 'none';
        }
    });

    // Update active state of pagination buttons
    const paginationButtons = container.querySelectorAll('.table-pagination button[data-page]');
    paginationButtons.forEach(button => {
        button.classList.toggle('active', parseInt(button.textContent) === pageNum);
    });

    // Reset the pagination based on the current page
    updatePagination(container, rows.length);
}

function resetTableSearch(container) {
    const table = container.querySelector('table');
    const rows = table.querySelectorAll('tbody tr');
    rows.forEach(row => {
        row.classList.remove('filtered-out');
        row.style.display = '';
    });
    updatePagination(container, rows.length);
    showPage(container, 1);
}

function initializeAllTableSearch() {
    document.querySelectorAll('.table-search').forEach(searchInput => {
        searchInput.addEventListener('input', function() {
            const container = this.closest('.card').querySelector('.table-container');
            filterTable(container, this.value);
        });
    });
}

// This should be called when the table is first displayed
function initializeTablePagination(container) {
    updatePagination(container);
    showPage(container, 1);
}

document.addEventListener('DOMContentLoaded', function() {
    console.log("DOM fully loaded");
    initializeAllTableSearch();
    // Other initialization code...
});


// Function to show a modal with an image
function showModal(imgSrc, button) {
    var modal = document.getElementById('figureModal');
    var modalImg = document.getElementById('modalImage');
    modal.style.display = "block";
    modalImg.src = imgSrc;

    // Add a class to the button when it's active
    button.classList.add('active-button');
}

// Function to close the modal
function closeModal() {
    var modal = document.getElementById('figureModal');
    modal.style.display = "none";

    // Remove the class from all buttons
    document.querySelectorAll('.toggle-button').forEach(button => {
        button.classList.remove('active-button');
    });
}

// Event listeners for modal behavior
document.addEventListener('DOMContentLoaded', function() {
    var modal = document.getElementById('figureModal');
    var span = document.getElementsByClassName("close")[0];

    // When the user clicks on <span> (x), close the modal
    span.onclick = closeModal;

    // When the user clicks anywhere outside of the modal content, close it
    window.onclick = function(event) {
        if (event.target == modal) {
            closeModal();
        }
    }

    // Ensure modal is hidden initially
    modal.style.display = "none";
});

// Function to initialize table pagination
function initializeTablePagination(tableContainer) {
    var tableContent = tableContainer.querySelector('.table-content');
    var table = tableContent.querySelector('table');
    var rows = table.querySelectorAll('tbody tr');
    var rowsPerPage = 10;
    var pageCount = Math.ceil(rows.length / rowsPerPage);
    var paginationContainer = tableContainer.querySelector('.table-pagination');
    var currentPage = 1;

    // Function to show a specific page of table rows
    function showPage(pageNum) {
        if (pageNum < 1 || pageNum > pageCount) return;

        currentPage = pageNum;
        var start = (pageNum - 1) * rowsPerPage;
        var end = start + rowsPerPage;

        rows.forEach((row, index) => {
            row.style.display = (index >= start && index < end) ? '' : 'none';
        });

        // Update the active class on the pagination buttons
        var paginationButtons = paginationContainer.querySelectorAll('button[data-page]');
        paginationButtons.forEach(button => {
            button.classList.remove('active');
        });
        paginationContainer.querySelector(`[data-page="${pageNum}"]`).classList.add('active');
    }

    // Function to create a pagination button
    function createPaginationButton(pageNum) {
        var btn = document.createElement('button');
        btn.innerText = pageNum;
        btn.setAttribute('data-page', pageNum);
        btn.addEventListener('click', function(e) {
            showPage(pageNum);
        });
        return btn;
    }

    // Function to setup pagination buttons for the table
    function setupPagination() {
        paginationContainer.innerHTML = '';

        // Previous page button
        var prevBtn = document.createElement('button');
        prevBtn.innerText = '<';
        prevBtn.addEventListener('click', function() {
            showPage(currentPage - 1);
        });
        paginationContainer.appendChild(prevBtn);

        // Page number buttons
        if (pageCount <= 5) {
            // If there are 5 or fewer pages, show all page numbers
            for (var i = 1; i <= pageCount; i++) {
                var btn = createPaginationButton(i);
                if (i === 1) {
                    btn.classList.add('active');
                }
                paginationContainer.appendChild(btn);
            }
        } else {
            // Show first three pages, ellipsis, and last page
            for (var i = 1; i <= 3; i++) {
                var btn = createPaginationButton(i);
                if (i === 1) {
                    btn.classList.add('active');
                }
                paginationContainer.appendChild(btn);
            }

            // Add ellipsis
            var ellipsis = document.createElement('span');
            ellipsis.innerText = '...';
            paginationContainer.appendChild(ellipsis);

            // Add last page
            var lastPageBtn = createPaginationButton(pageCount);
            paginationContainer.appendChild(lastPageBtn);
        }

        // Next page button
        var nextBtn = document.createElement('button');
        nextBtn.innerText = '>';
        nextBtn.addEventListener('click', function() {
            showPage(currentPage + 1);
        });
        paginationContainer.appendChild(nextBtn);

        // Page number input
        var pageInput = document.createElement('input');
        pageInput.type = 'number';
        pageInput.min = 1;
        pageInput.max = pageCount;
        pageInput.placeholder = 'Go ';
        pageInput.classList.add('pagination-input');
        pageInput.addEventListener('keypress', function(e) {
            if (e.key === 'Enter') {
                var pageNum = parseInt(pageInput.value);
                if (!isNaN(pageNum) && pageNum >= 1 && pageNum <= pageCount) {
                    showPage(pageNum);
                }
            }
        });
        paginationContainer.appendChild(pageInput);
    }

    setupPagination();
    showPage(1); // Show the first page by default
}


// Function to initialize toggle-button 
document.querySelectorAll('.toggle-button').forEach(button => {
    button.addEventListener('click', function() {
        // Check if the button already has the active class
        if (this.classList.contains('active-button')) {
            // If it does, remove the active class to "close" it
            this.classList.remove('active-button');
        } else {
            // Remove active class from all buttons
            document.querySelectorAll('.toggle-button').forEach(btn => btn.classList.remove('active-button'));

            // Add active class to the clicked button to "open" it
            this.classList.add('active-button');
        }
    });
});
