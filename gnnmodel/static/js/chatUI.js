// For the "ongoing_turn" action handler
function showGeneratingMessage() {
  // Create container for the "generating" message
  const generatingContainer = document.getElementById("bottom-chat-log");

  // Limpa o conteúdo anterior
  generatingContainer.innerHTML = "";

  // Create message bubble
  const messageBubble = document.createElement("div");
  messageBubble.className = "p-3 ms-3 bot-message d-flex align-items-center";

  // Bootstrap spinner
  const spinner = document.createElement("div");
  spinner.className = "spinner-border spinner-border-sm text-primary me-2";
  spinner.setAttribute("role", "status");
  spinner.innerHTML = '<span class="visually-hidden">Loading...</span>';

  // Create message text
  const messageText = document.createElement("p");
  messageText.className = "small mb-0 me-3";
  messageText.textContent = "Generating response...";

  // Botão de parar
  const stopButton = document.createElement("button");
  stopButton.className = "btn btn-sm btn-danger";
  stopButton.textContent = "Stop";
  stopButton.onclick = function () {
    chatSocket.send(JSON.stringify({ action: "stop_generating" }));
    stopButton.disabled = true;
    stopButton.textContent = "Stopping...";
  };

  // Monta os elementos
  messageBubble.appendChild(spinner);
  messageBubble.appendChild(messageText);
  messageBubble.appendChild(stopButton);
  generatingContainer.appendChild(messageBubble);

  // Scroll to the bottom
  generatingContainer.scrollIntoView();
}

function updateActiveModel(modelName) {
  // Remove active class from all items
  var modelItems = document.querySelectorAll("#models-list .dropdown-item");
  modelItems.forEach(function (item) {
    item.classList.remove("active");
  });

  // Add active class to the selected model
  var modelItems = document.querySelectorAll("#models-list .dropdown-item");
  modelItems.forEach(function (item) {
    if (item.textContent === modelName) {
      item.classList.add("active");
    }
  });
}

// Add a function to populate the models dropdown
function populateModelsList(models, currentModel) {
  var modelsList = document.getElementById("models-list");
  modelsList.innerHTML = "";

  models.forEach(function (model) {
    var li = document.createElement("li");
    var a = document.createElement("a");
    a.className = "dropdown-item" + (model === currentModel ? " active" : "");
    a.href = "#";
    a.textContent = model;
    a.onclick = function () {
      changeModel(model);
      return false;
    };
    li.appendChild(a);
    modelsList.appendChild(li);
  });
}

// Update the chat log display
function updateChatLog() {
  const chatLog = document.querySelector("#chat-log");

  if (!chatLog) return; // Ensure chatLog exists

  // Clear existing content
  chatLog.innerHTML = "";

  // Create DOM elements for each message
  messages.forEach(function (msg) {
    // Create main container
    const messageContainer = document.createElement("div");
    messageContainer.className = "d-flex flex-row mb-4";

    // Set alignment based on message source
    if (msg.source === "assistant") {
      messageContainer.classList.add("justify-content-start");
    } else {
      messageContainer.classList.add("justify-content-end");
    }

    // Create message bubble
    const messageBubble = document.createElement("div");
    messageBubble.className = "p-3 ms-3";

    // Set bubble style based on message source
    if (msg.source === "assistant") {
      messageBubble.classList.add("bot-message");
    } else {
      messageBubble.classList.add("user-message");
    }

    // Create message text
    const messageText = document.createElement("p");
    messageText.className = "small mb-0";
    // Check if there's file info to display
    if (msg.file_info) {
      if (msg.file_info.type.startsWith("image/")) {
        // Display image preview
        messageText.innerHTML = `<img src="${msg.file_info.data}" alt="${msg.file_info.name}" style="max-width: 200px; max-height: 200px; display: block; margin-bottom: 5px;">`;
      } else if (msg.file_info.type === "application/pdf") {
        // Display PDF icon/placeholder
        messageText.innerHTML = `<i class="fas fa-file-pdf me-2"></i><span>${msg.file_info.name}</span><br>`;
      } else {
        // Display generic file info
        messageText.innerHTML = `<i class="fas fa-file me-2"></i><span>${msg.file_info.name}</span><br>`;
      }
    }
    messageText.innerHTML += msg.msg; // Append the text message

    // Assemble the elements
    messageBubble.appendChild(messageText);
    messageContainer.appendChild(messageBubble);
    chatLog.appendChild(messageContainer);
  });

  // Create container for the "generating" message
  const generatingContainer = document.createElement("div");
  generatingContainer.className = "d-flex flex-row mb-4 justify-content-center";
  generatingContainer.id = "bottom-chat-log";
  chatLog.appendChild(generatingContainer);
  generatingContainer.scrollIntoView();
}

// Populate the sessions dropdown
function populateSessionsList(sessions) {
  var sessionsList = document.getElementById("sessions-list");
  sessionsList.innerHTML = "";

  if (sessions.length === 0) {
    var li = document.createElement("li");
    li.innerHTML =
      '<span class="dropdown-item text-muted">No saved sessions</span>';
    sessionsList.appendChild(li);
  } else {
    sessions.forEach(function (session) {
      var li = document.createElement("li");
      li.className = "d-flex align-items-center";

      // Create session link
      var a = document.createElement("a");
      a.className = "dropdown-item flex-grow-1";
      a.href = "#";
      a.textContent = session.name;
      a.dataset.sessionId = session.session_id;
      a.onclick = function () {
        loadSession(session.session_id, session.name);
        return false;
      };

      // Create delete button
      var deleteBtn = document.createElement("button");
      deleteBtn.type = "button";
      deleteBtn.className = "btn btn-sm text-danger me-2";
      deleteBtn.innerHTML = '<i class="fas fa-trash"></i>';
      deleteBtn.title = "Delete session";
      deleteBtn.onclick = function (e) {
        e.stopPropagation(); // Prevent dropdown item click
        showDeleteConfirmation(session.session_id, session.name);
        return false;
      };

      li.appendChild(a);
      li.appendChild(deleteBtn);
      sessionsList.appendChild(li);
    });
  }

  // Add divider and option to create new session
  var divider = document.createElement("li");
  divider.innerHTML = '<hr class="dropdown-divider">';
  sessionsList.appendChild(divider);

  var newSessionLi = document.createElement("li");
  var newSessionA = document.createElement("a");
  newSessionA.className = "dropdown-item text-success";
  newSessionA.href = "#";
  newSessionA.textContent = "Create New Session";
  newSessionA.onclick = function () {
    showNewSessionModal();
    return false;
  };
  newSessionLi.appendChild(newSessionA);
  sessionsList.appendChild(newSessionLi);
}

function populateToolsList(tools) {
  var toolsList = document.getElementById("tools-list");
  toolsList.innerHTML = "";

  var li = document.createElement("li");
  var confirmBtn = document.createElement("button");
  confirmBtn.className = "btn btn-sm btn-outline-success rounded-circle mx-1";
  confirmBtn.innerHTML = "<i class='fas fa-check'></i>";
  confirmBtn.title = "Confirm";
  confirmBtn.type = "button";
  confirmBtn.onclick = function (event) {
    chatSocket.send(
      JSON.stringify({
        action: "change_tools",
        tools: selectedTools,
      })
    );
  };
  li.appendChild(confirmBtn);

  var selectAllBtn = document.createElement("button");
  selectAllBtn.className =
    "btn btn-sm btn-outline-secondary rounded-circle mx-1";
  selectAllBtn.innerHTML = "<i class='fas fa-list-check'></i>";
  selectAllBtn.title = "Select/Deselect all tools";
  selectAllBtn.type = "button";
  selectAllBtn.onclick = function (event) {
    event.stopPropagation();
    var checkboxes = toolsList.querySelectorAll("input[type='checkbox']");
    var allChecked = Array.from(checkboxes).every(
      (checkbox) => checkbox.checked
    );
    checkboxes.forEach((checkbox) => {
      checkbox.checked = !allChecked;
      if (checkbox.checked) {
        if (!selectedTools.includes(checkbox.value)) {
          selectedTools.push(checkbox.value);
        }
      } else {
        selectedTools = selectedTools.filter((t) => t !== checkbox.value);
      }
    });
  };
  li.appendChild(selectAllBtn);
  toolsList.appendChild(li);

  // Add divider
  var divider = document.createElement("li");
  divider.innerHTML = '<hr class="dropdown-divider">';
  toolsList.appendChild(divider);

  tools.forEach(function (tool) {
    var li = document.createElement("li");
    li.onclick = function (event) {
      event.stopPropagation();
    };
    li.className = "d-flex align-items-center";

    var label = document.createElement("label");
    label.className = "dropdown-item";

    var checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.value = tool;
    checkbox.checked = selectedTools.includes(tool);
    checkbox.name = "tool-item";
    checkbox.onchange = function (event) {
      if (this.checked) {
        if (!selectedTools.includes(tool)) {
          selectedTools.push(tool);
        }
      } else {
        selectedTools = selectedTools.filter((t) => t !== tool);
      }
    };

    label.appendChild(checkbox);
    label.appendChild(document.createTextNode(" " + tool));

    toolInfoBtn = document.createElement("button");
    toolInfoBtn.className = "btn btn-sm btn-outline-info rounded-circle mx-1";
    toolInfoBtn.innerHTML = "<i class='fas fa-info-circle'></i>";
    toolInfoBtn.title = "Info";
    toolInfoBtn.type = "button";
    toolInfoBtn.setAttribute("data-bs-toggle", "modal");
    toolInfoBtn.setAttribute("data-bs-target", "#toolInfoModal");
    toolInfoBtn.onclick = function (event) {
      event.stopPropagation();
      showToolInfo(tool);
    };
    li.appendChild(label);
    li.appendChild(toolInfoBtn);
    toolsList.appendChild(li);
  });
}

function showToolInfo(toolName) {
  let modal = document.getElementById("toolInfoModal");
  modal.onclick = (event) => {
    event.stopPropagation();
  };

  // Define o conteúdo do modal
  document.getElementById(
    "tool-info-body"
  ).innerHTML = `<strong>${toolName}</strong><br>${
    toolDescriptions[toolName] || "No description available."
  }`;
}

// Simple toast notification function
function showToast(message, type = "success") {
  // Create toast container if it doesn't exist
  let toastContainer = document.getElementById("toast-container");
  if (!toastContainer) {
    toastContainer = document.createElement("div");
    toastContainer.id = "toast-container";
    toastContainer.className = "position-fixed top-0 end-0 p-3";
    document.body.appendChild(toastContainer);
  }

  // Create toast
  const toastId = "toast-" + Date.now();
  const toast = document.createElement("div");
  toast.id = toastId;
  toast.className = `toast align-items-center ${
    type === "error" ? "bg-danger" : "bg-success"
  } text-white`;
  toast.setAttribute("role", "alert");
  toast.setAttribute("aria-live", "assertive");
  toast.setAttribute("aria-atomic", "true");

  toast.innerHTML = `
    <div class="d-flex">
      <div class="toast-body">
        ${message}
      </div>
      <button type="button" class="btn-close btn-close-white me-2 m-auto" data-bs-dismiss="toast" aria-label="Close"></button>
    </div>
  `;

  toastContainer.appendChild(toast);

  // Show toast
  const bsToast = new bootstrap.Toast(toast, {
    autohide: true,
    delay: 5000,
  });
  bsToast.show();

  // Remove toast after it's hidden
  toast.addEventListener("hidden.bs.toast", function () {
    toast.remove();
  });
}

function removeSelectedFile() {
  selectedFile = null;
  document.getElementById("file-input").value = ""; // Clear the file input
  const previewContainer = document.getElementById("file-preview-container");
  previewContainer.style.display = "none";
  previewContainer.innerHTML = "";
}
