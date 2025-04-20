var wss_protocol = window.location.protocol == "https:" ? "wss://" : "ws://";
var chatSocket;
var messages = [];
var currentSessionId = null;
var currentSessionName = "New Session";
var availableModels = [];
var currentModelName = "";

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
      console.log("Requesting sessions list");
      chatSocket.send(
        JSON.stringify({
          action: "get_sessions",
        })
      );
    }, 500); // Small delay to ensure connection is fully established
  };

  chatSocket.onmessage = function (e) {
    var data = JSON.parse(e.data);
    console.log("Received message:", data);

    // Handle different types of messages
    if (data.action) {
      console.log("Handling action:", data.action);
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
var deleteSessionId = null;
var deleteSessionName = null;

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
      break;
    case "model_changed":
      // Update the current model when changed
      if (data.model_name) {
        currentModelName = data.model_name;
        document.getElementById("current-model-name").textContent =
          currentModelName;
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
      document.getElementById("bottom-chat-log").innerHTML = "";
      document.getElementById("bottom-chat-log").scrollIntoView();
      break;
    case "ongoing_turn":
      chat_log_end = ` <div class="d-flex flex-row mb-4 justify-content-center" id="bottom-chat-log">
  <div class="p-3 ms-3 bot-message">
      <p class="small mb-0">Generating response...</p></div></div>`;
      document.querySelector("#chat-log").innerHTML += chat_log_end;
      document.getElementById("bottom-chat-log").scrollIntoView();
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
  var str = "";
  messages.forEach(function (msg) {
    str += `<div class="d-flex flex-row mb-4 ${
      msg.source == "assistant"
        ? "justify-content-start"
        : "justify-content-end"
    }">
          <div class="p-3 ms-3 ${
            msg.source == "assistant" ? "bot-message" : "user-message"
          }">
          <p class="small mb-0">${msg.msg}</p></div></div>`;
  });

  document.querySelector("#chat-log").innerHTML = str;
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
  document.querySelector("#chat-message-input").focus();
  document.querySelector("#chat-message-input").onkeyup = function (e) {
    if (e.keyCode === 13) {
      // enter, return
      document.querySelector("#chat-message-submit").click();
    }
  };

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
