var selectAllItems = "#select-all";
var checkboxItem = ":checkbox";

$(selectAllItems).click(function() {

  if (this.checked) {
    $(checkboxItem).each(function() {
      this.checked = true;
    });
  } else {
    $(checkboxItem).each(function() {
      this.checked = false;
    });
  }

});