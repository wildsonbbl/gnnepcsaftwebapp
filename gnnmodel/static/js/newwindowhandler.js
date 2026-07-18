// Armazena a função original do window.open
const originalWindowOpen = window.open;

// Função auxiliar para verificar de forma robusta se um link é externo
function isLinkExternal(url) {
  try {
    // Resolve URLs relativas (como "/pure/") para absolutas com base no endereço atual
    const absoluteUrl = new URL(url, window.location.href);

    // Se não for HTTP ou HTTPS (ex: mailto:, tel:), tratamos como externo/especial
    if (absoluteUrl.protocol !== "http:" && absoluteUrl.protocol !== "https:") {
      return true;
    }

    // Compara os domínios/hostnames (ignora a porta para evitar problemas com portas dinâmicas)
    const targetHost = absoluteUrl.hostname.replace("www.", "");
    const currentHost = window.location.hostname.replace("www.", "");

    return targetHost !== currentHost;
  } catch (e) {
    return true;
  }
}

// Função principal que decide como abrir o link
function handleOpenLink(url) {
  // Resolve a URL para garantir que passamos o caminho completo ao Python
  const absoluteUrl = new URL(url, window.location.href).href;

  if (
    window.pywebview &&
    window.pywebview.api &&
    window.pywebview.api.open_link
  ) {
    console.log("Direcionando link ao pywebview: ", absoluteUrl);
    window.pywebview.api.open_link(absoluteUrl);
  } else {
    console.log("pywebview API não está pronta, usando fallback");
    originalWindowOpen(absoluteUrl, "_blank");
  }
}

// Sobrescreve o window.open global do navegador
window.open = function (url, target, features) {
  console.log("Link aberto via window.open: ", url);
  handleOpenLink(url);
  return null;
};

// Intercepta cliques globais em links <a>
document.addEventListener("click", function (e) {
  const target = e.target.closest("a");
  if (target) {
    const href = target.getAttribute("href");
    if (!href || href.startsWith("#") || href.startsWith("javascript:")) {
      return; // Ignora âncoras internas da página ou scripts
    }

    const absoluteUrl = new URL(href, window.location.href).href;
    const isExternal = isLinkExternal(absoluteUrl);
    const isModifierClick = e.ctrlKey || e.metaKey || e.shiftKey;
    const isOpenInNewTab = target.target === "_blank" || isModifierClick;

    // Intercepta APENAS se for link externo OU se foi explicitamente configurado para abrir em nova aba
    if (isExternal || isOpenInNewTab) {
      e.preventDefault();
      console.log("Clique interceptado (externo ou _blank): ", absoluteUrl);
      handleOpenLink(absoluteUrl);
    } else {
      // Links internos normais (como a navegação do menu) continuam o fluxo normal do navegador
      console.log("Navegação interna permitida: ", absoluteUrl);
    }
  }
});

// Log para confirmar que a API do pywebview foi devidamente carregada
window.addEventListener("pywebviewready", function () {
  console.log("Integração do pywebview ativada com sucesso!");
});

// Função auxiliar para ler o valor atual de um cookie
function getCookie(name) {
  let cookieValue = null;
  if (document.cookie && document.cookie !== "") {
    const cookies = document.cookie.split(";");
    for (let i = 0; i < cookies.length; i++) {
      const cookie = cookies[i].trim();
      if (cookie.substring(0, name.length + 1) === name + "=") {
        cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
        break;
      }
    }
  }
  return cookieValue;
}

// Intercepta o envio de formulários e atualiza o token CSRF com o valor mais recente do Cookie
document.addEventListener("submit", function (e) {
  const csrfInput = e.target.querySelector('input[name="csrfmiddlewaretoken"]');
  if (csrfInput) {
    const currentCookieToken = getCookie("csrftoken");
    if (currentCookieToken && csrfInput.value !== currentCookieToken) {
      console.log(
        "Sincronizando token CSRF desatualizado com o cookie mais recente...",
      );
      csrfInput.value = currentCookieToken;
    }
  }
});
