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
var selectedFile = null;
var mcpConfigPath = "";
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

// Event listeners
document.addEventListener("DOMContentLoaded", function () {
  // Initialize chat
  initializeChat();

  // Set up input and send button

  const textarea = document.getElementById("chat-message-input");
  const fileInput = document.getElementById("file-input"); // Get file input
  const attachButton = document.getElementById("attach-file-button"); // Get attach button
  const filePreviewContainer = document.getElementById(
    "file-preview-container"
  );
  filePreviewContainer.style.display = "none"; // Hide initially

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

  // Trigger file input click when attach button is clicked
  attachButton.onclick = function () {
    fileInput.click();
  };

  // Handle file selection
  fileInput.onchange = function (e) {
    if (e.target.files.length > 0) {
      const file = e.target.files[0];
      const maxSize = 50 * 1024 * 1024; // Example: 50MB limit

      if (file.size > maxSize) {
        showToast(
          `File is too large. Maximum size is ${maxSize / 1024 / 1024} MB.`,
          "error"
        );
        removeSelectedFile();
        return;
      }

      const reader = new FileReader();

      reader.onload = function (event) {
        selectedFile = {
          name: file.name,
          type: file.type,
          data: event.target.result, // Base64 data URL
        };
        // Show preview
        filePreviewContainer.innerHTML = `
                <span class="me-2">${file.name} (${(file.size / 1024).toFixed(
          1
        )} KB)</span>
                <button type="button" class="btn-close btn-sm" aria-label="Remove file" onclick="removeSelectedFile()"></button>
            `;
        filePreviewContainer.style.display = "block";
      };
      reader.onerror = function (event) {
        console.error("File could not be read! Error: ", event.target.error);
        showToast("Error reading file.", "error");
        removeSelectedFile();
      };
      reader.readAsDataURL(file); // Read as Base64 Data URL
    } else {
      removeSelectedFile();
    }
  };

  // Add listener for the MCP Config button
  const mcpConfigBtn = document.getElementById("config-mcp-btn");
  if (mcpConfigBtn) {
    mcpConfigBtn.addEventListener("click", function () {
      // Request current MCP config content when modal is about to be shown
      chatSocket.send(JSON.stringify({ action: "get_mcp_config_content" }));
      const mcpConfigError = document.getElementById("mcp-config-error");
      mcpConfigError.classList.add("d-none"); // Clear previous errors
      mcpConfigError.textContent = "";
    });
  } else {
    console.error("MCP Config button not found.");
  }

  // Add listener for the Update Ollama Models button
  const updateOllamaModelsBtn = document.getElementById(
    "update-ollama-models-btn"
  );
  if (updateOllamaModelsBtn) {
    updateOllamaModelsBtn.onclick = function () {
      showToast("Updating Ollama models...", "info");
      chatSocket.send(JSON.stringify({ action: "update_ollama_models" }));
    };
  }

  // Add listener for the Save MCP Config button in the modal
  const saveMcpConfigBtn = document.getElementById("save-mcp-config-btn");
  if (saveMcpConfigBtn) {
    saveMcpConfigBtn.onclick = function () {
      const newConfigContent = document.getElementById(
        "mcp-config-content-input"
      ).value;
      const mcpConfigError = document.getElementById("mcp-config-error");
      mcpConfigError.classList.add("d-none");
      mcpConfigError.textContent = "";

      try {
        JSON.parse(newConfigContent); // Basic validation for JSON
      } catch (e) {
        mcpConfigError.textContent = "Invalid JSON format: " + e.message;
        mcpConfigError.classList.remove("d-none");
        showToast("MCP Configuration is not valid JSON.", "error");
        return;
      }

      chatSocket.send(
        JSON.stringify({
          action: "save_mcp_config_content",
          content: newConfigContent,
        })
      );
    };
  } else {
    console.error("Save MCP Config button not found.");
  }

  document.querySelector("#chat-message-submit").onclick = function (e) {
    var messageInputDom = document.querySelector("#chat-message-input");
    var message = messageInputDom.value;

    // Check if there's either text or a file selected
    if (!message.trim() && !selectedFile) return;

    const messagePayload = {
      text: message,
    };

    if (selectedFile) {
      messagePayload.file = selectedFile;
    }

    chatSocket.send(JSON.stringify(messagePayload));

    messageInputDom.value = "";
    removeSelectedFile(); // Clear file after sending
    messageInputDom.style.height = "auto"; // Reset textarea height
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
