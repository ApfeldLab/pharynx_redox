TITLES = getList("image.titles");

Dialog.create("Create Masks");
Dialog.addMessage("Select the image to create mask.");
Dialog.addMessage("Rotations will be based on this mask.");
Dialog.addChoice("Image", TITLES);
Dialog.show()