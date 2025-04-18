var wss_protocol = window.location.protocol == "https:" ? "wss://" : "ws://";
var chatSocket;
var messages = [];
var currentSessionId = null;
var currentSessionName = "New Session";

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
  }
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
      var a = document.createElement("a");
      a.className = "dropdown-item";
      a.href = "#";
      a.textContent = session.name;
      a.dataset.sessionId = session.session_id;
      a.onclick = function () {
        loadSession(session.session_id, session.name);
        return false;
      };
      li.appendChild(a);
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
});
