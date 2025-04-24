var wss_protocol = window.location.protocol == "https:" ? "wss://" : "ws://";
var chatSocket;
var messages = [];
var currentSessionId = null;
var currentSessionName = "New Session";
var availableModels = [];
var currentModelName = "";
var deleteSessionId = null;
var deleteSessionName = null;
var availableTools = [];
var selectedTools = [];
const toolDescriptions = {
  ToolA: "Descrição da ToolA.",
  ToolB: "Descrição da ToolB.",
  ToolC: "Descrição da ToolC.",
  // Adicione mais ferramentas aqui
};

// Initialize the chat
function initializeChat(sessionId = null) {
  // Close existing socket if any
  if (chatSocket) {
    chatSocket.close();
  }

  // Create WebSocket connection with session ID if provided
  var wsUrl = wss_protocol + window.location.host + "/ws/chat/";
  if (sessionId) {
    wsUrl += sessionId + "/";
  }

  chatSocket = new WebSocket(wsUrl);
  setupChatSocketHandlers();

  // Update UI with current session name
  document.getElementById("current-session-name").textContent =
    currentSessionName;
}

// Set up WebSocket event handlers
function setupChatSocketHandlers() {
  chatSocket.onopen = function (e) {
    console.log("Connected to chat server");
    // Request list of sessions
    setTimeout(function () {
      chatSocket.send(
        JSON.stringify({
          action: "get_sessions",
        })
      );
    }, 500); // Small delay to ensure connection is fully established
  };

  chatSocket.onmessage = function (e) {
    var data = JSON.parse(e.data);

    // Handle different types of messages
    if (data.action) {
      handleActionMessage(data);
    } else if (data.text) {
      handleChatMessage(data.text);
    }
  };

  chatSocket.onclose = function (e) {
    console.log("Socket closed unexpectedly, please reload the page.");
  };
}

// Handle action messages (session management)
function handleActionMessage(data) {
  switch (data.action) {
    case "sessions_list":
      populateSessionsList(data.sessions);
      break;
    case "session_loaded":
      // Set the current session when initially loaded
      currentSessionId = data.session_id;
      currentSessionName = data.name;
      document.getElementById("current-session-name").textContent =
        currentSessionName;

      // Handle model information
      if (data.model_name) {
        currentModelName = data.model_name;
        document.getElementById("current-model-name").textContent =
          currentModelName;
      }
      // Populate available models if provided
      if (data.available_models && Array.isArray(data.available_models)) {
        availableModels = data.available_models;
        populateModelsList(availableModels, currentModelName);
      }
      if (data.selected_tools && Array.isArray(data.selected_tools)) {
        selectedTools = [...data.selected_tools];
      } else if (data.available_tools && Array.isArray(data.available_tools)) {
        selectedTools = [...data.available_tools];
      }
      // Populate available tools if provided
      if (data.available_tools && Array.isArray(data.available_tools)) {
        availableTools = data.available_tools;
        populateToolsList(availableTools);
      }
      if (data.tool_descriptions) {
        Object.assign(toolDescriptions, data.tool_descriptions);
      }
      break;
    case "tools_changed":
      if (data.selected_tools) {
        selectedTools = [...data.selected_tools];
        populateToolsList(availableTools);
        showToast("Tools changed successfully");
      }
      break;
    case "model_changed":
      // Update the current model when changed
      if (data.model_name) {
        currentModelName = data.model_name;
        document.getElementById("current-model-name").textContent =
          currentModelName;
        // Update the active class in the models dropdown
        updateActiveModel(currentModelName);
      }
      break;
    case "session_created":
      currentSessionId = data.session_id;
      currentSessionName = data.name;
      document.getElementById("current-session-name").textContent =
        currentSessionName;
      // Clear messages for new session
      messages = [];
      updateChatLog();
      // Close the modal
      var modal = bootstrap.Modal.getInstance(
        document.getElementById("newSessionModal")
      );
      if (modal) modal.hide();
      // Atualiza a lista de sessões
      chatSocket.send(
        JSON.stringify({
          action: "get_sessions",
        })
      );
      break;
    case "load_messages":
      messages = data.messages;
      updateChatLog();
      break;
    case "session_renamed":
      if (data.session_id === currentSessionId) {
        currentSessionName = data.name;
        document.getElementById("current-session-name").textContent =
          currentSessionName;
      }
      // Refresh sessions list
      chatSocket.send(
        JSON.stringify({
          action: "get_sessions",
        })
      );
      // Close the modal
      var modal = bootstrap.Modal.getInstance(
        document.getElementById("renameSessionModal")
      );
      if (modal) modal.hide();
      break;
    case "end_turn":
      const generatingContainer = document.getElementById("bottom-chat-log");
      generatingContainer.innerHTML = "";
      generatingContainer.scrollIntoView();
      if (data.type === "stop_action") {
        showToast("Stopped generating", "error");
      }
      if (data.type === "interrupted") {
        showToast("LLM interrupted generating", "error");
      }
      break;
    case "ongoing_turn":
      showGeneratingMessage();
      break;
    case "session_deleted":
      if (data.success) {
        // If the current session was deleted, update UI
        if (data.session_id === currentSessionId) {
          messages = [];
          updateChatLog();
        }

        // Close the modal
        var modal = bootstrap.Modal.getInstance(
          document.getElementById("deleteSessionModal")
        );
        if (modal) modal.hide();

        // Show success message
        showToast("Session deleted successfully");
      } else {
        showToast("Failed to delete session", "error");
      }
      break;
  }
}

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

// Add a function to change the model
function changeModel(modelName) {
  if (modelName === currentModelName) return;

  chatSocket.send(
    JSON.stringify({
      action: "change_model",
      model_name: modelName,
      tools: selectedTools,
    })
  );
}

// Handle chat messages
function handleChatMessage(message) {
  if ((message.source == "assistant") | (message.source == "user")) {
    messages.push(message);
  }
  updateChatLog();
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
    messageText.innerHTML = msg.msg;

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

// Function to show delete confirmation modal
function showDeleteConfirmation(sessionId, sessionName) {
  deleteSessionId = sessionId;
  deleteSessionName = sessionName;

  document.getElementById("delete-session-name").textContent = sessionName;
  var modal = new bootstrap.Modal(
    document.getElementById("deleteSessionModal")
  );
  modal.show();
}

// Function to delete a session
function deleteSession() {
  if (!deleteSessionId) return;

  chatSocket.send(
    JSON.stringify({
      action: "delete_session",
      session_id: deleteSessionId,
    })
  );
}

// Simple toast notification function
function showToast(message, type = "success") {
  // Create toast container if it doesn't exist
  let toastContainer = document.getElementById("toast-container");
  if (!toastContainer) {
    toastContainer = document.createElement("div");
    toastContainer.id = "toast-container";
    toastContainer.className = "position-fixed bottom-0 end-0 p-3";
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
    delay: 3000,
  });
  bsToast.show();

  // Remove toast after it's hidden
  toast.addEventListener("hidden.bs.toast", function () {
    toast.remove();
  });
}

// Load a session
function loadSession(sessionId, sessionName) {
  currentSessionId = sessionId;
  currentSessionName = sessionName;
  document.getElementById("current-session-name").textContent =
    currentSessionName;

  chatSocket.send(
    JSON.stringify({
      action: "load_session",
      session_id: sessionId,
    })
  );
}

// Show the new session modal
function showNewSessionModal() {
  var modal = new bootstrap.Modal(document.getElementById("newSessionModal"));
  document.getElementById("new-session-name").value = "";
  modal.show();
}

// Show the rename session modal
function showRenameSessionModal() {
  if (!currentSessionId) {
    alert("Please create or select a session first");
    return;
  }

  var modal = new bootstrap.Modal(
    document.getElementById("renameSessionModal")
  );
  document.getElementById("rename-session-name").value = currentSessionName;
  modal.show();
}

// Create a new session
function createNewSession() {
  var sessionName =
    document.getElementById("new-session-name").value.trim() || "New Session";

  chatSocket.send(
    JSON.stringify({
      action: "create_session",
      name: sessionName,
      model_name: currentModelName,
      tools: selectedTools,
    })
  );
}

// Rename the current session
function renameCurrentSession() {
  if (!currentSessionId) return;

  var newName = document.getElementById("rename-session-name").value.trim();
  if (!newName) return;

  chatSocket.send(
    JSON.stringify({
      action: "rename_session",
      session_id: currentSessionId,
      name: newName,
    })
  );
}

// Event listeners
document.addEventListener("DOMContentLoaded", function () {
  // Initialize chat
  initializeChat();

  // Set up input and send button

  const textarea = document.getElementById("chat-message-input");
  textarea.focus();
  textarea.onkeyup = function (e) {
    if (e.key === "Enter" && !e.shiftKey) {
      // enter, return
      document.querySelector("#chat-message-submit").click();
      this.style.height = "auto";
    }
  };

  textarea.addEventListener("input", function () {
    this.style.height = "auto";
    this.style.height = Math.min(this.scrollHeight, 300) + "px"; // 300px é o limite máximo
  });

  document.querySelector("#chat-message-submit").onclick = function (e) {
    var messageInputDom = document.querySelector("#chat-message-input");
    var message = messageInputDom.value;
    if (!message.trim()) return;

    chatSocket.send(
      JSON.stringify({
        text: message,
      })
    );

    messageInputDom.value = "";
  };

  // Export chat log
  document
    .getElementById("chat-log-save_button")
    .addEventListener("click", function () {
      const dataStr =
        "data:text/json;charset=utf-8," +
        encodeURIComponent(JSON.stringify(messages));
      const downloadAnchorNode = document.createElement("a");
      downloadAnchorNode.setAttribute("href", dataStr);
      downloadAnchorNode.setAttribute(
        "download",
        `${currentSessionName || "chat"}.json`
      );
      document.body.appendChild(downloadAnchorNode);
      downloadAnchorNode.click();
      downloadAnchorNode.remove();
    });

  // New session button
  document
    .getElementById("new-session-btn")
    .addEventListener("click", showNewSessionModal);

  // Create session button in modal
  document
    .getElementById("create-session-btn")
    .addEventListener("click", createNewSession);

  // Rename session button
  document
    .getElementById("rename-session-btn")
    .addEventListener("click", showRenameSessionModal);

  // Save rename button in modal
  document
    .getElementById("save-rename-btn")
    .addEventListener("click", renameCurrentSession);
  // Add event listener for delete confirmation
  document
    .getElementById("confirm-delete-btn")
    .addEventListener("click", deleteSession);
});
