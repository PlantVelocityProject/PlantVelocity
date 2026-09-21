// Search page titles in the sidebar; page content can be added independently.
const search = document.getElementById('navigation-search');
const navigation = document.getElementById('documentation-nav');
if (search && navigation) {
  search.hidden = false;
  search.addEventListener('input', () => {
    const query = search.value.trim().toLowerCase();
    let visibleCount = 0;
    navigation.querySelectorAll('.nav-group').forEach(group => {
      let matches = 0;
      group.querySelectorAll('a').forEach(link => {
        link.hidden = !link.textContent.toLowerCase().includes(query);
        if (!link.hidden) matches++;
      });
      group.hidden = matches === 0;
      visibleCount += matches;
    });
    navigation.querySelector('.search-status').hidden = visibleCount !== 0;
  });
}
