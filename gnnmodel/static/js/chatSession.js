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
